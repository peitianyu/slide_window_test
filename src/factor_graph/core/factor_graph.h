#ifndef __FACTOR_GRAPH_H__
#define __FACTOR_GRAPH_H__

#include "factor.h"
#include "variable.h"
#include "pattern_ordering.h"
#include "core/tt_assert.h"

#include <map>
#include <vector>
#include <Eigen/Dense>
struct VariablePattern
{
    int total_variables_dim = 0;
    std::map<Variable *, PatternOrdering> variable_lookup_table;
    Eigen::MatrixXd H;
    Eigen::VectorXd b;

    inline const PatternOrdering &VariableLookup(Variable *var){
        tt_assert(variable_lookup_table.count(var) != 0 && "Variable not found in variable_lookup_table");
        return variable_lookup_table[var];
    }
};

class FactorGraph
{
public:
    struct Option
    {
        int max_iterations;        // Maximum number of iterations.
        double relative_error_th; // Maximum relative error decrease.
        double absolute_error_th; // Maximum absolute error decrease.

        Option(int max_iter = 20, double rel_err_th = 1e-5, double abs_err_th = 1e-5)
            : max_iterations(max_iter), relative_error_th(rel_err_th), absolute_error_th(abs_err_th) {}
    };

    FactorGraph(Option option = Option()) : m_option(option) {}
    
    ~FactorGraph()
    {
        for (Factor* f : m_factors){delete f;}
        for (Variable* v : m_variables){delete v;}
    }
    
    void AddFactor(Factor* f) { m_factors.push_back(f); }
    void AddVariable(Variable* v) { m_variables.push_back(v); }

    const std::vector<Variable*>& GetVariables() const{ return m_variables; }
    const std::vector<Factor*>& GetFactors() const{return m_factors; }

    bool OptimizeGN(std::vector<Factor*> factors)
    {
        for(Factor* f : factors) ConstructSparsityPattern(f);

        for(auto& var : m_optimized_variables){
            var->Print();
        }
           
        double current_error = 0.5 * ComputeErrorNormSquared(factors); // 用于收敛判断
        double last_error = current_error;

        bool converged = false;
        for(int iter = 0; iter < m_option.max_iterations; ++iter)
        {
            Eigen::VectorXd dx = m_variable_pattern.H.ldlt().solve(m_variable_pattern.b);
            UpdateSparsityPattern(dx);
            std::cout << "GN Iteration: " << iter << "current_error: " << current_error << std::endl;
            current_error = 0.5 * ComputeErrorNormSquared(factors);
            if(IsConverged(iter, last_error, current_error)){
                converged = true;
                break;
            }
            last_error = current_error;
        }

        std::cout << "-------------------------" << std::endl;
        std::cout << "converged: " << converged << std::endl;
        for(auto& var : m_optimized_variables)
            var->Print();

        return converged;

        // std::cout << "H: \n" << m_variable_pattern.H << std::endl;
        // std::cout << "b: \n" << m_variable_pattern.b.transpose() << std::endl;

        // Variable* marg_var = factors.front()->VariableAt(0);
        // MarginalizeSparsityPattern(marg_var);

        // std::cout << "H_marg: \n" << m_variable_pattern.H << std::endl;
        // std::cout << "b_marg: \n" << m_variable_pattern.b.transpose() << std::endl;
    }
private:
    double ComputeErrorNormSquared(const std::vector<Factor*>& factors)
    {
        double error_norm_squared = 0.0;
        for(Factor* f : factors)
            error_norm_squared += f->Error().squaredNorm();
        return error_norm_squared;
    }

    bool IsConverged(int iter_num, double current_error, double new_error)
    {
        if (!std::isfinite(new_error)) return false;

        if (iter_num == m_option.max_iterations) return false;

        const double error_decrease = current_error - new_error;
        if (error_decrease > m_option.absolute_error_th) return false;
        if((error_decrease / current_error) > m_option.relative_error_th) return false;
           
        return true;
    }

    void ConstructSparsityPattern(Factor* f)
    {
        // 1. 构造PatternOrdering
        int total_variables_dim = m_variable_pattern.total_variables_dim;
        for(int i = 0; i < f->NumVariables(); ++i)
        {
            Variable* var = f->VariableAt(i);
            if(m_variable_pattern.variable_lookup_table.count(var) == 0){
                m_optimized_variables.push_back(var);
                m_variable_pattern.variable_lookup_table[var] = PatternOrdering().SetDim(var->Dim()).SetIdx(m_variable_pattern.total_variables_dim);
                total_variables_dim += var->Dim();
            }
        }

        if(total_variables_dim > m_variable_pattern.total_variables_dim){
            Eigen::VectorXd b = Eigen::VectorXd::Zero(total_variables_dim);
            b.head(m_variable_pattern.total_variables_dim) = m_variable_pattern.b;
            m_variable_pattern.b = b;
            
            Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_variables_dim, total_variables_dim);
            H.topLeftCorner(m_variable_pattern.total_variables_dim, m_variable_pattern.total_variables_dim) = m_variable_pattern.H;
            m_variable_pattern.H = H;

            m_variable_pattern.total_variables_dim = total_variables_dim;
        }

        // 2. 构造H和b
        const int n_rows = f->ErrorDim();
        const int num_vars = f->NumVariables();
        std::vector<int> vars_cols(num_vars, -1);
        std::vector<int> vars_dim(num_vars, -1);
        int n_cols = 0;
        for (int i = 0; i < num_vars; ++i)
        {
            Variable *var = f->VariableAt(i);
            if (!var->fixed){
                const int var_dim = m_variable_pattern.VariableLookup(var).Dim();
                vars_cols[i] = n_cols;
                vars_dim[i] = var_dim;
                n_cols += var_dim;
            }
        }

        // 2.1 构造H
        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(n_rows, n_cols);
        for (int i = 0, start_col = 0; i < num_vars; ++i)
        {
            if (!f->VariableAt(i)->fixed)
            {
                Eigen::MatrixXd J = f->Jacobian(i);
                Js.block(0, start_col, J.rows(), J.cols()) = J;
                start_col += J.cols();
            }
        }
        const Eigen::MatrixXd JtJ = Js.transpose() * Js;
        for(int i = 0; i < num_vars; ++i){
            if(f->VariableAt(i)->fixed) continue;
            const int var_idx = m_variable_pattern.variable_lookup_table[f->VariableAt(i)].Idx();
            const int var_dim = vars_dim[i];
            m_variable_pattern.H.block(var_idx, var_idx, var_dim, var_dim) += JtJ.block(vars_cols[i], vars_cols[i], var_dim, var_dim);
        }

        // 2.2 构造b
        const Eigen::VectorXd Jtb = Js.transpose() * f->Error();
        for(int i = 0; i < num_vars; ++i){
            if(f->VariableAt(i)->fixed) continue;
            const int var_idx = m_variable_pattern.variable_lookup_table[f->VariableAt(i)].Idx();
            const int var_dim = vars_dim[i];
            m_variable_pattern.b.segment(var_idx, var_dim) -= Jtb.segment(vars_cols[i], var_dim);
        }
    }

    void MarginalizeSparsityPattern(Variable *marg_var)
    {
        // 1. 重新构造variable_lookup_table, 将marg_var放到最后
        int marg_var_idx = m_variable_pattern.variable_lookup_table[marg_var].Idx();
        int marg_var_dim = m_variable_pattern.variable_lookup_table[marg_var].Dim();
        for(auto& it : m_variable_pattern.variable_lookup_table){
            if(it.second.Idx() > marg_var_idx){
                it.second.SetIdx(it.second.Idx() - marg_var_dim);
            }
        }
        m_variable_pattern.variable_lookup_table[marg_var].SetIdx(m_variable_pattern.total_variables_dim - marg_var_dim);

        // 2. 重新构造H和b, 将marg_var对应的行和列放到最后

        int idx = marg_var_idx;  // marg_var的索引
        int dim = marg_var_dim;  // marg_var的维度
        int reserve_size = m_variable_pattern.total_variables_dim;  // 矩阵的维度

        // 将 row i 移动矩阵最下边
        Eigen::MatrixXd temp_rows = m_variable_pattern.H.block(idx, 0, dim, reserve_size);
        Eigen::MatrixXd temp_botRows = m_variable_pattern.H.block(idx + dim, 0, reserve_size - idx - dim, reserve_size);
        m_variable_pattern.H.block(idx, 0, reserve_size - idx - dim, reserve_size) = temp_botRows;
        m_variable_pattern.H.block(reserve_size - dim, 0, dim, reserve_size) = temp_rows;

        // 将 col i 移动矩阵最右边
        Eigen::MatrixXd temp_cols = m_variable_pattern.H.block(0, idx, reserve_size, dim);
        Eigen::MatrixXd temp_rightCols = m_variable_pattern.H.block(0, idx + dim, reserve_size, reserve_size - idx - dim);
        m_variable_pattern.H.block(0, idx, reserve_size, reserve_size - idx - dim) = temp_rightCols;
        m_variable_pattern.H.block(0, reserve_size - dim, reserve_size, dim) = temp_cols;

        // 3. 重新构造H和b
        Eigen::MatrixXd H11 = m_variable_pattern.H.block(0, 0, reserve_size - dim, reserve_size - dim);
        Eigen::MatrixXd H12 = m_variable_pattern.H.block(0, reserve_size - dim, reserve_size - dim, dim);
        Eigen::MatrixXd H21 = m_variable_pattern.H.block(reserve_size - dim, 0, dim, reserve_size - dim);
        Eigen::MatrixXd H22 = m_variable_pattern.H.block(reserve_size - dim, reserve_size - dim, dim, dim);
        Eigen::VectorXd b1 = m_variable_pattern.b.segment(0, reserve_size - dim);
        Eigen::VectorXd b2 = m_variable_pattern.b.segment(reserve_size - dim, dim);

        /*
            // 简化运算, 所以并不直接使用
            Eigen::MatrixXd H22_inv = H22.inverse();
            Eigen::MatrixXd H_schur = H11 - H12 * H22_inv * H21;
            Eigen::VectorXd b_schur = b1 - H12 * H22_inv * b2;
            m_variable_pattern.H = H_schur;
            m_variable_pattern.b = b_schur;
        */
        Eigen::MatrixXd temp = H12 * H22.inverse();
        m_variable_pattern.H = H11 - temp * H21;
        m_variable_pattern.b = b1 - temp * b2;
        
        // 4. 更新total_variables_dim
        m_variable_pattern.total_variables_dim -= marg_var_dim;

        // 5. 更新variable_lookup_table
        m_variable_pattern.variable_lookup_table.erase(marg_var);
    }

    void UpdateSparsityPattern(const Eigen::VectorXd& dx)
    {
        for (int i = 0, d = 0, count = m_optimized_variables.size(); i < count; ++i){
            Variable *var = m_optimized_variables[i];
            if (!var->fixed){
                const int vd = m_variable_pattern.variable_lookup_table[var].Dim();
                var->Plus(dx.segment(d, vd));
                d += vd;
            }
        }
    }
private:
    Option m_option;
    std::vector<Variable*> m_variables;
    std::vector<Variable*> m_optimized_variables;
    std::vector<Factor*> m_factors;
    VariablePattern m_variable_pattern;  
};


#endif // __FACTOR_GRAPH_H__