#ifndef __FACTOR_GRAPH_H__
#define __FACTOR_GRAPH_H__

#include "factor.h"
#include "variable.h"
#include "pattern_ordering.h"

#include <map>
#include <vector>
#include <Eigen/Sparse>

struct SparsityPattern
{
    std::map<Variable *, PatternOrdering> variable_lookup_table;
    Eigen::MatrixXd H;
    Eigen::VectorXd b;

    inline const PatternOrdering &VariableLookup(Variable *var){
        tt_assert(variable_lookup_table.count(var) != 0);
        return variable_lookup_table[var];
    }
};

class FactorGraph
{
public:
    FactorGraph() = default;
    
    ~FactorGraph()
    {
        for (Factor* f : m_factors){delete f;}
        for (Variable* v : m_variables){delete v;}
    }
    
    void AddFactor(Factor* f) { m_factors.push_back(f); }
    void AddVariable(Variable* v) { m_variables.push_back(v); }

    void AddFactors(const std::vector<Factor*>& factors) { m_factors.insert(m_factors.end(), factors.begin(), factors.end()); }
    void AddVariables(const std::vector<Variable*>& variables) { m_variables.insert(m_variables.end(), variables.begin(), variables.end()); }

    const std::vector<Variable*>& GetVariables() const{ return m_variables; }
    const std::vector<Factor*>& GetFactors() const{return m_factors; }

    void OptimizeGN()
    {

    }
private:
    void ConstructSparsityPattern(Factor* f)
    {

    }

    void MarginalizeSparsityPattern(Variable *marg_var)
    {

    }

    void UpdateSparsityPattern(const Eigen::VectorXd& dx)
    {

    }
private:
    std::vector<Variable*> m_variables;
    std::vector<Factor*> m_factors;
    SparsityPattern m_sparsity_pattern;
};

#endif // __FACTOR_GRAPH_H__