#ifndef __SLIDE_WINDOW_OPTIMIZE_H__
#define __SLIDE_WINDOW_OPTIMIZE_H__

#include <iostream>
#include <map>
#include "factor_graph.h"
#include <Eigen/Sparse>
#include "utils.h"
#include "sparsity_pattern.h"

class SlideWindowOptimize
{
public:
    struct Option
    {
        int slide_window_size; // Maximum number of iterations.
        int max_iterations;        // Maximum number of iterations.
        double relative_error_th; // Maximum relative error decrease.
        double absolute_error_th; // Maximum absolute error decrease.

        Option(int slide_window_size = 1000, int max_iterations = 3000, double relative_error_th = 1e-5, double absolute_error_th = 1e-5)
            : slide_window_size(slide_window_size), max_iterations(max_iterations), relative_error_th(relative_error_th), absolute_error_th(absolute_error_th)
        {}
    };

    SlideWindowOptimize(Option option) : m_option(option) {}

    void AddFactor(Factor *factor) { m_graph.AddFactor(factor);}
    void AddVariable(Variable *variable) { m_graph.AddVariable(variable);}

    void AddFactors(const std::vector<Factor*>& factors) {m_graph.AddFactors(factors);}
    void AddVariables(const std::vector<Variable*>& variables) {m_graph.AddVariables(variables);}

    bool OptimizeGN()
    {
        GRAPH_LOG("Started optimization with %i factors and %i variables.", (int)graph->GetFactors().size(), (int)graph->GetVariables().size());
        if(m_graph.GetVariables().size() < m_option.slide_window_size)
            m_sparsity_pattern_builder.ConstructSparsityPattern(m_graph, &m_pattern);  // 初始化稀疏图缓存器

        double current_error = 0.5 * ComputeErrorNormSquared(m_graph); // 用于收敛判断
        GRAPH_LOG("Initial error: %f", current_error);
        double last_error = current_error;

        bool converged = false;
        for(int iter = 0; iter < m_option.max_iterations; ++iter)
        {
            if (!Iterate(&m_graph, &m_pattern))
                return false;

            current_error = 0.5 * ComputeErrorNormSquared(m_graph);
            GRAPH_LOG("Iteration %i: error %f", iter, current_error);

            if (!ContinueIteratingCheck(iter, last_error, current_error, &converged))
                break;
            
            last_error = current_error;
        }

        if(m_graph.GetVariables().size() >= m_option.slide_window_size)
            m_sparsity_pattern_builder.Marginalize(m_graph, &m_pattern, m_graph.GetVariables().front());

        return converged;
    }
private:
    bool Iterate(FactorGraph *graph, SparsityPattern *pattern);

    void LinearizeSingleFactor(Factor *factor, SparsityPattern *pattern);

    void Marginalize();
private:
    Option m_option;
    SparsityPatternBuilder  m_sparsity_pattern_builder;
    SparsityPattern m_pattern;
    FactorGraph m_graph;
};





#endif // __SLIDE_WINDOW_OPTIMIZE_H__