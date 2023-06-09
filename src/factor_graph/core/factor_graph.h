#ifndef __FACTOR_GRAPH_H__
#define __FACTOR_GRAPH_H__

#include "factor.h"
#include "variable.h"
#include "pattern_ordering.h"
#include "core/tt_assert.h"

#include <map>
#include <vector>
#include <Eigen/Sparse>

struct SparsityPattern
{
    int total_variables_dim = 0;
    std::map<Variable *, PatternOrdering> variable_lookup_table;
    std::vector<Eigen::Triplet<double>> triplet_list; // H稀疏矩阵的三元组表示, 仅存储下三角部分
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
        // 1. 构造PatternOrdering
        int total_variables_dim = m_sparsity_pattern.total_variables_dim;
        for(int i = 0; i < f->NumVariables(); ++i)
        {
            Variable* var = f->VariableAt(i);
            if(m_sparsity_pattern.variable_lookup_table.count(var) == 0){
                m_sparsity_pattern.variable_lookup_table[var] = PatternOrdering().SetDim(var->Dim()).SetIdx(m_sparsity_pattern.total_variables_dim);
                total_variables_dim += var->Dim();
            }
        }

        if(total_variables_dim > m_sparsity_pattern.total_variables_dim){
            Eigen::VectorXd b = Eigen::VectorXd::Zero(total_variables_dim);
            m_sparsity_pattern.b.head(m_sparsity_pattern.total_variables_dim) = b;
            m_sparsity_pattern.total_variables_dim = total_variables_dim;
        }
        

        // 2. 构造triplet_list和b
        const int n_rows = f->ErrorDim();
        const int num_vars = f->NumVariables();
        std::vector<int> vars_cols(num_vars, -1);
        std::vector<int> vars_dim(num_vars, -1);
        int n_cols = 0;
        for (int i = 0; i < num_vars; ++i)
        {
            Variable *var = f->VariableAt(i);
            if (!var->fixed){
                const int var_dim = m_sparsity_pattern.VariableLookup(var).Dim();
                vars_cols[i] = n_cols;
                vars_dim[i] = var_dim;
                n_cols += var_dim;
            }
        }
        
        tt_assert(n_cols != 0, "All variables connected to this factor are fixed.");

        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(n_rows, n_cols);
        for (int i = 0, start_col = 0; i < num_vars; ++i) {
            if (!f->VariableAt(i)->fixed) {
                Eigen::MatrixXd J = f->Jacobian(i);
                Js.block(0, start_col, J.rows(), J.cols()) = J;
                start_col += J.cols();
            }
        }

        // 求解JtJ获得整个H,之后会对稀疏H矩阵进行赋值
        Eigen::MatrixXd JtJ = Js.transpose() * Js;

        // 2.1 构造triplet_list
        // 只填充下三角
        for(int i = 0; i < num_vars; ++i)
        {
            if(f->VariableAt(i)->fixed) continue;

            const int H_col = m_sparsity_pattern.VariableLookup(f->VariableAt(i)).Idx();
            const int JtJ_col = vars_cols[i];
            for(int j = i; j < num_vars; ++j)
            {
                if(f->VariableAt(j)->fixed) continue;

                const int H_row = m_sparsity_pattern.VariableLookup(f->VariableAt(j)).Idx();
                const int JtJ_row = vars_cols[j];
                for (int JtJ_i = JtJ_col, H_i = H_col; JtJ_i < (JtJ_col + vars_dim[i]); ++JtJ_i, ++H_i){
                    for (int JtJ_j = JtJ_row, H_j = H_row; JtJ_j < (JtJ_row + vars_dim[j]); ++JtJ_j, ++H_j){
                        m_sparsity_pattern.triplet_list.push_back(Eigen::Triplet<double>(H_i, H_j, JtJ(JtJ_i, JtJ_j)));
                    }
                } 
            }     
        }

        // 2.2 构造b
        const Eigen::VectorXd Jtb = Js.transpose() * f->Error();
        for(int i = 0; i < num_vars; ++i){
            if(f->VariableAt(i)->fixed) continue;
            const int var_idx = m_sparsity_pattern.VariableLookup(f->VariableAt(i)).Idx();
            const int var_dim = vars_dim[i];
            m_sparsity_pattern.b.segment(var_idx, var_dim) -= Jtb.segment(vars_cols[i], var_dim);
        }
    }

    void MarginalizeSparsityPattern(Variable *marg_var)
    {
        // 1. 重新构造variable_lookup_table, 将marg_var放到最后
        int marg_var_idx = m_sparsity_pattern.variable_lookup_table[marg_var].Idx();
        int marg_var_dim = m_sparsity_pattern.variable_lookup_table[marg_var].Dim();
        for(auto& it : m_sparsity_pattern.variable_lookup_table){
            if(it.second.Idx() > marg_var_idx){
                it.second.SetIdx(it.second.Idx() - marg_var_dim);
            }
        }
        m_sparsity_pattern.variable_lookup_table[marg_var].SetIdx(m_sparsity_pattern.total_variables_dim - marg_var_dim);

        // 2. 重新构造triplet_list和b, 将marg_var对应的行和列放到最后
        Eigen::SparseMatrix<double> H_marg = Eigen::SparseMatrix<double>(m_sparsity_pattern.total_variables_dim, m_sparsity_pattern.total_variables_dim);
        H_marg.setFromTriplets(m_sparsity_pattern.triplet_list.begin(), m_sparsity_pattern.triplet_list.end());

        int idx = marg_var_idx;  // marg_var的索引
        int dim = marg_var_dim;  // marg_var的维度
        int reserve_size = m_sparsity_pattern.total_variables_dim;  // 矩阵的维度

        // 将 row i 移动矩阵最下边
        Eigen::MatrixXd temp_rows = H_marg.block(idx, 0, dim, reserve_size);
        Eigen::MatrixXd temp_botRows = H_marg.block(idx + dim, 0, reserve_size - idx - dim, reserve_size);
        H_marg.block(idx, 0, reserve_size - idx - dim, reserve_size) = temp_botRows;
        H_marg.block(reserve_size - dim, 0, dim, reserve_size) = temp_rows;

        // 将 col i 移动矩阵最右边
        Eigen::MatrixXd temp_cols = H_marg.block(0, idx, reserve_size, dim);
        Eigen::MatrixXd temp_rightCols = H_marg.block(0, idx + dim, reserve_size, reserve_size - idx - dim);
        H_marg.block(0, idx, reserve_size, reserve_size - idx - dim) = temp_rightCols;
        H_marg.block(0, reserve_size - dim, reserve_size, dim) = temp_cols;

        // 3. 重新构造triplet_list和b
        Eigen::MatrixXd H11 = H_marg.block(0, 0, reserve_size - dim, reserve_size - dim);
        Eigen::MatrixXd H12 = H_marg.block(0, reserve_size - dim, reserve_size - dim, dim);
        Eigen::MatrixXd H21 = H_marg.block(reserve_size - dim, 0, dim, reserve_size - dim);
        Eigen::MatrixXd H22 = H_marg.block(reserve_size - dim, reserve_size - dim, dim, dim);
        Eigen::VectorXd b1 = m_sparsity_pattern.b.head(reserve_size - dim);

        Eigen::MatrixXd H22_inv = H22.inverse();

        Eigen::MatrixXd H_schur = H11 - H12 * H22_inv * H21;
        Eigen::VectorXd b_schur = b1 - H12 * H22_inv * b1;

        // 重新构造triplet_list, 仅保留下三角
        m_sparsity_pattern.triplet_list.clear();
        for(int i = 0; i < H_schur.rows(); ++i){
            for(int j = 0; j < H_schur.cols(); ++j){
                m_sparsity_pattern.triplet_list.push_back(Eigen::Triplet<double>(i, j, H_schur(i, j)));
            }
        }

        m_sparsity_pattern.b = b_schur;
        m_sparsity_pattern.total_variables_dim -= marg_var_dim;

        // 4. 重新构造variable_lookup_table, 删除marg_var
        m_sparsity_pattern.variable_lookup_table.erase(marg_var);
    }

    void UpdateSparsityPattern(const Eigen::VectorXd& dx)
    {
        for (int i = 0, d = 0, count = m_variables.size(); i < count; ++i){
            Variable *var = m_variables[i];
            if (!var->fixed){
                const int vd = m_sparsity_pattern.variable_lookup_table[var].Dim();
                var->Plus(dx.segment(d, vd));
                d += vd;
            }
        }
    }
private:
    std::vector<Variable*> m_variables;
    std::vector<Factor*> m_factors;
    SparsityPattern m_sparsity_pattern;
};

#endif // __FACTOR_GRAPH_H__