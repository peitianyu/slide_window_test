#ifndef __FACTOR_H__
#define __FACTOR_H__ 

#include <vector>
#include <array>
#include <memory>
#include <Eigen/Core>
#include "variable.h"
#include "core/tt_assert.h"
#include <iostream>

class Factor
{
public:
    static constexpr int kMaxVariables = 2;

    Factor() = default;

    virtual ~Factor() {}

    int NumVariables() const { return m_num_variables; }

    void AddVariable(Variable *v){
        tt_assert(m_num_variables < kMaxVariables && "Too many variables in this factor.");
        m_variables[m_num_variables] = v;
        ++m_num_variables;
    }

    Variable *VariableAt(int idx) const {
        tt_assert(m_num_variables > 0 && idx < m_num_variables && "Invalid variable index.");
        return m_variables[idx];
    }

    // Information matrix. Defaults to identity.
    void SetInfoMatrix(const Eigen::MatrixXd &info_matrix) { m_info_matrix = info_matrix; }
    Eigen::MatrixXd GetSqrtInfoMatrix() const { return m_info_matrix; }

    // Dimensionality of the error.
    virtual int ErrorDim() const = 0;
    virtual Eigen::VectorXd Error() const = 0;

    // Returns e1 - e2.
    virtual Eigen::VectorXd SubtractError(const Eigen::VectorXd &e1, const Eigen::VectorXd &e2) const = 0;

    // Jacobian wrt to the variable at idx. 
    // Defaults to computing the jacobian numerically.
    virtual Eigen::MatrixXd Jacobian(int idx) const {
        tt_assert(m_num_variables > 0 && idx < m_num_variables && "Invalid variable index.");
        return ComputeNumericalJacobian(m_variables[idx]);
    }
protected:
    Eigen::MatrixXd ComputeNumericalJacobian(Variable * v) const
    {
        constexpr double h = 1e-5;
        const int N = v->Dim();
        const int M = this->ErrorDim();
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(M, N);
        Eigen::VectorXd dx = Eigen::VectorXd::Zero(N);
        Eigen::VectorXd dy0 = this->Error();
        constexpr double k = 1.0 / (2.0 * h);
        for (int i = 0; i < N; ++i)
        {
            dx(i) = h;
            v->Plus(dx); // right
            const Eigen::VectorXd dy1 = this->SubtractError(this->Error(), dy0);
            dx(i) = -2.0 * h;
            v->Plus(dx); // left
            const Eigen::VectorXd dy2 = this->SubtractError(this->Error(), dy0);
            dx(i) = h;
            v->Plus(dx); // return to original state.
            dx(i) = 0.0;
            J.col(i) << (dy1 - dy2) * k;
        }
        return J;
    }
private:
    std::array<Variable *, kMaxVariables> m_variables;
    int m_num_variables = 0;
    Eigen::MatrixXd m_info_matrix;
};

#endif // __FACTOR_H__