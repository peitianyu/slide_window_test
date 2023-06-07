#include"graph_optimize.h"

GraphOptimize::GraphOptimize(Option option)
    : m_option(option)
{}

bool GraphOptimize::OptimizeGN(FactorGraph *graph)
{
    // 1. 构造稀疏图缓存器
    GRAPH_LOG("Started optimization with %i factors and %i variables.", (int)graph->GetFactors().size(), (int)graph->GetVariables().size());
    SparsityPattern pattern;
    m_sparsity_pattern_builder.ConstructSparsityPattern(*graph, &pattern);  // 初始化稀疏图缓存器

    double current_error = 0.5 * ComputeErrorNormSquared(*graph); // 用于收敛判断
    GRAPH_LOG("Initial error: %f", current_error);
    double last_error = current_error;

    bool converged = false;
    for(int iter = 0; iter < m_option.max_iterations; ++iter)
    {
        if (!Iterate(graph, &pattern))
            return false;

        current_error = 0.5 * ComputeErrorNormSquared(*graph);
        GRAPH_LOG("Iteration %i: error %f", iter, current_error);

        if (!ContinueIteratingCheck(iter, last_error, current_error, &converged))
            break;
        
        last_error = current_error;
    }

    return converged;
}

bool GraphOptimize::Iterate(FactorGraph *graph, SparsityPattern *pattern)
{
    const std::vector<Factor *> &factors = graph->GetFactors();
    pattern->b.setZero();
    pattern->triplet_list.clear();
    for(Factor *factor: factors)
        LinearizeSingleFactor(factor, pattern);// 结算H与b, 并赋值给pattern
    GRAPH_ASSERT(!pattern->triplet_list.empty());
    pattern->H.setFromTriplets(pattern->triplet_list.begin(), pattern->triplet_list.end());// 从triplet_list获取H矩阵

    // TODO: 如果H并不是稀疏矩阵,通过近似最小度排序变为稀疏矩阵
    // AMD 或者 COLAMD
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(pattern->H);
    Eigen::VectorXd dx = chol.solve(pattern->b);
    if (chol.info() != Eigen::Success)
    {
        GRAPH_LOG("Cholesky decomposition failed.");
        return false;
    }
    m_sparsity_pattern_builder.UpdateSparsityPattern(graph, pattern, dx);
    return true;
}

void GraphOptimize::LinearizeSingleFactor(Factor *factor, SparsityPattern *pattern)
{
    const int n_rows = factor->ErrorDim();
    const int num_variables = factor->NumVariables();
    std::vector<int> vars_cols(num_variables, -1);
    std::vector<int> vars_dim(num_variables, -1);
    int n_cols = 0;
    for (int i = 0; i < num_variables; ++i)
    {
        Variable *var = factor->VariableAt(i);
        if (!var->fixed)
        {
            const int var_dim = pattern->VariableLookup(var).Dim();
            vars_cols[i] = n_cols;
            vars_dim[i] = var_dim;
            n_cols += var_dim;
        }
        
    }

    if (n_cols == 0)
    {
        GRAPH_LOG("All variables connected to this factor are fixed.");
        return;
    }
    Eigen::MatrixXd Js(n_rows, n_cols);
    for (int i = 0, start_col = 0; i < num_variables; ++i)
    {
        if (!factor->VariableAt(i)->fixed)
        {
            Eigen::MatrixXd J = factor->Jacobian(i);
            Js.block(0, start_col, J.rows(), J.cols()) = J;
            start_col += J.cols();
        }
    }

    // 求解JtJ获得整个H,之后会对稀疏H矩阵进行赋值
    Eigen::MatrixXd JtJ(Js.cols(), Js.cols());
    JtJ.noalias() = Js.transpose() *  factor->GetSqrtInfoMatrix() * Js; // Eigen默认会解决混淆问题，如果你确定不会出现混淆，可以使用noalias（）来提效

    // 现在我们需要将贡献值添加到h中，注意我们只填充下三角形部分。 
    // 通过雅可比求解H矩阵
    for (int i = 0; i < num_variables; ++i)
    {
        if (factor->VariableAt(i)->fixed)
            continue;

        const int H_col = pattern->VariableLookup(factor->VariableAt(i)).Idx();
        const int JtJ_col = vars_cols[i];
        for (int j = i; j < num_variables; ++j)
        {
            if (factor->VariableAt(j)->fixed)
                continue;
            
            const int H_row = pattern->VariableLookup(factor->VariableAt(j)).Idx();
            const int JtJ_row = vars_cols[j];
            for (int JtJ_i = JtJ_col, H_i = H_col; JtJ_i < (JtJ_col + vars_dim[i]); ++JtJ_i, ++H_i){
                for (int JtJ_j = JtJ_row, H_j = H_row; JtJ_j < (JtJ_row + vars_dim[j]); ++JtJ_j, ++H_j){
                    pattern->triplet_list.push_back(SparsityPattern::T(H_j, H_i, JtJ(JtJ_j, JtJ_i)));
                }
            }
        }
    }
    
    // Handle the vector b.
    // 对b赋值
    const Eigen::VectorXd Jtb = Js.transpose() *  factor->GetSqrtInfoMatrix() * factor->Error();
    for (int i = 0; i < num_variables; ++i)
    {
        Variable *var = factor->VariableAt(i);
        if (!var->fixed)
        {
            const int var_idx = pattern->VariableLookup(var).Idx();
            const int var_dim = vars_dim[i];
            pattern->b.segment(var_idx, var_dim) -= Jtb.segment(vars_cols[i], var_dim);
        }
    }
}

bool GraphOptimize::ContinueIteratingCheck(int iter_num, double current_error, double new_error, bool *converged)
{
    if (!std::isfinite(new_error)){
        return false;
    }

    if (iter_num == m_option.max_iterations){
        GRAPH_LOG("Max iterations reached.");
        return false;
    }

    const double error_decrease = current_error - new_error;
    if ((error_decrease <= m_option.absolute_error_th) || ((error_decrease / current_error) <= m_option.relative_error_th))
    {
        GRAPH_LOG("Converged.");
        *converged = true;
        return false;
    }

    return true;
}

double GraphOptimize::ComputeErrorNormSquared(const FactorGraph &graph)
{
    const std::vector<Factor *> &factors = graph.GetFactors();
    double error = 0.0;
    for(Factor * factor : factors)
        error += factor->Error().squaredNorm();
    return error;
}
