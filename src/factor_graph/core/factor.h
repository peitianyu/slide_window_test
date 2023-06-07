#ifndef __FACTOR_H__
#define __FACTOR_H__

#include<vector>
#include<memory>
#include<Eigen/Core>
#include"variable.h"
#include "utils.h"
#include <iostream>

class Factor
{
public:
    static constexpr int kMaxVariables = 2;

    Factor() = default;

    int NumVariables() const { return m_num_variables; }

    void AddVariable(Variable *v);

    Variable *VariableAt(int idx) const
    {
        GRAPH_ASSERT(m_num_variables > 0 && idx < m_num_variables);
        return m_variables[idx];
    }

    // Dimensionality of the error.
    virtual int ErrorDim() const = 0;

    virtual Eigen::VectorXd Error() const = 0;

    // Returns e1 - e2.
    virtual Eigen::VectorXd SubtractError(const Eigen::VectorXd &e1, const Eigen::VectorXd &e2) const = 0;

    // Jacobian wrt to the variable at idx. 
    // Defaults to computing the jacobian numerically.
    virtual Eigen::MatrixXd Jacobian(int idx) const
    {
        GRAPH_ASSERT(m_num_variables > 0 && idx < m_num_variables);
        return ComputeNumericalJacobian(m_variables[idx]);
    }
    
    void SetInfoMatrix(const Eigen::MatrixXd &info_matrix) { m_info_matrix = info_matrix; }

    Eigen::MatrixXd GetSqrtInfoMatrix() const { return m_info_matrix; }
protected:
    Eigen::MatrixXd ComputeNumericalJacobian(Variable * v) const;
    std::array<Variable *, kMaxVariables> m_variables;
    int m_num_variables = 0;
    Eigen::MatrixXd m_info_matrix;
};

#endif // __FACTOR_H__