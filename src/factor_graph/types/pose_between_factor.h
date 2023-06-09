#ifndef __POSE_BETWEEN_FACTOR_H__
#define __POSE_BETWEEN_FACTOR_H__

#include "../core/factor.h"
#include "pose2d.h"

class PoseBetweenFactor : public Factor
{
public:
    PoseBetweenFactor(Pose2d *v_a, Pose2d *v_b, double x_ab, double y_ab,
                 double yaw_ab_rad, const Eigen::Matrix3d &info_matrix) : m_pos_ab(x_ab, y_ab), m_yaw_ab_rad(yaw_ab_rad)
    {
        AddVariable(v_a);
        AddVariable(v_b);
        SetInfoMatrix(info_matrix);
    }

    virtual int ErrorDim() const { return 3; }

    virtual Eigen::VectorXd Error() const override
    {
        tt_assert(this->NumVariables() == 2);
        const Pose2d *v_a = static_cast<Pose2d *>(this->VariableAt(0));
        const Pose2d *v_b = static_cast<Pose2d *>(this->VariableAt(1));
        Eigen::Vector3d r;
        Eigen::Vector2d pos_ab_pred = {v_b->x() - v_a->x(), v_b->y() - v_a->y()};
        r.head<2>() = Eigen::Rotation2Dd(v_a->yaw_rad()).toRotationMatrix().transpose() * pos_ab_pred - m_pos_ab;
        r(2) = NormalizeAngle((v_b->yaw_rad() - v_a->yaw_rad()) - m_yaw_ab_rad);
        return r;
    }

    virtual Eigen::VectorXd SubtractError(const Eigen::VectorXd &e1, const Eigen::VectorXd &e2) const override
    {
        Eigen::Vector3d diff;
        diff << (e1(0) - e2(0)), (e1(1) - e2(1)), NormalizeAngle(e1(2) - e2(2));
        return diff;
    }
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
private:
    double NormalizeAngle(double theta_rad) const
    {
        // Normalize the angle to the range [-pi, pi).
        constexpr double kPI = 3.14159265358979323846;
        constexpr double k2PI = 2.0 * kPI;
        return (theta_rad - k2PI * std::floor((theta_rad + kPI) / k2PI));
    }
private:
    Eigen::Vector2d m_pos_ab;
    double m_yaw_ab_rad;
};


#endif // __POSE_BETWEEN_FACTOR_H__