#include <ceres/ceres.h>

#include "my_types.h"
#include "utilities.hpp"

class AutoDiffFunctor {
public:
    AutoDiffFunctor(const myEdge& edge): measurement(edge.measurement), sqrt_info_matrix(edge.infoMatrix.array().sqrt()) {}

    template <typename T>
    bool operator()(const T* const v1, const T* const v2, T* residual) const {
        
        
        Eigen::Matrix<T, 3, 1> v1_mat{ v1[0], v1[1], v1[2] };
        Eigen::Matrix<T, 3, 1> v2_mat{ v2[0], v2[1], v2[2] };
        Eigen::Matrix<T, 3, 1> m = measurement.template cast<T>();
        Eigen::Map<Eigen::Matrix<T, 3, 1>> error{ residual };

        // use sqrt(info) * Z^(-1) * X1^(-1) & X2, involving two inverse operation
        // Eigen::Matrix<T, 3, 3> X1 = PoseToTrans(v1_mat);
        // Eigen::Matrix<T, 3, 3> X2 = PoseToTrans(v2_mat);
        // Eigen::Matrix<T, 3, 3> Z = PoseToTrans(m);

        // Eigen::Matrix<T, 3, 3> E = sqrt_info_matrix.template cast<T>() * Z.inverse() * (X1.inverse() * X2);
        // error = TransToPose(E);

        // calculate error from translation and rotation respectively and combine them together
        Eigen::Matrix<T, 3, 3> X1 = PoseToTrans(v1_mat);
        Eigen::Matrix<T, 3, 3> X2 = PoseToTrans(v2_mat);
        Eigen::Matrix<T, 3, 3> Z = PoseToTrans(m);
        
        Eigen::Matrix<T, 2, 2> Ri = X1.block(0, 0, 2, 2);
        Eigen::Matrix<T, 2, 2> Rj = X2.block(0, 0, 2, 2);
        Eigen::Matrix<T, 2, 2> Rij = Z.block(0, 0, 2, 2);

        Eigen::Matrix<T, 2, 1> ti{ v1_mat(0), v1_mat(1) };
        Eigen::Matrix<T, 2, 1> tj{ v2_mat(0), v2_mat(1) };
        Eigen::Matrix<T, 2, 1> tij{ m(0), m(1) };

        Eigen::Matrix<T, 2, 2> dRiT_dtheta;       //  derivative of Ri^T over theta
        dRiT_dtheta(0, 0) = T(-1) * Ri(1, 0); //  cosX -> -sinX
        dRiT_dtheta(0, 1) = T( 1) * Ri(0, 0); //  sinX ->  cosX
        dRiT_dtheta(1, 0) = T(-1) * Ri(0, 0); // -sinX -> -cosX
        dRiT_dtheta(1, 1) = T(-1) * Ri(1, 0); //  cosX -> -sinX

        // calcuate error & normalize error on theta
        error.template segment<2>(0) = Rij.transpose() * (Ri.transpose() * (tj - ti) - tij);
        error(2) = v2_mat(2) - v1_mat(2) - m(2);
        if (error(2) > T(M_PI)) {
            error(2) -= T(2 * M_PI);
        } else if (error(2) < T(-1 * M_PI)) {
            error(2) += T(2 * M_PI);
        }

        error = sqrt_info_matrix.template cast<T>() * error;

        // Eigen::Matrix<T, 3, 3> E = sqrt_info_matrix.template cast<T>() * Z.inverse() * (X1.inverse() * X2);
        // error = TransToPose(E);

        return true;
    }

    static ceres::CostFunction* create(const myEdge& edge) {
        return (new ceres::AutoDiffCostFunction<AutoDiffFunctor, 3, 3, 3>(new AutoDiffFunctor(edge)));
    }

private:
    Eigen::Vector3d measurement;
    Eigen::Matrix3d sqrt_info_matrix;
};

class AnalyticDiffFunction : public ceres::SizedCostFunction<3, 3, 3> {
public:
    virtual ~AnalyticDiffFunction() {}
    
    AnalyticDiffFunction(const myEdge& edge): measurement(edge.measurement), sqrt_info_matrix(edge.infoMatrix.array().sqrt()) {}
    
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {

        Eigen::Vector3d xi{ parameters[0][0], parameters[0][1], parameters[0][2] };
        Eigen::Vector3d xj{ parameters[1][0], parameters[1][1], parameters[1][2] };
        Eigen::Map<Eigen::Vector3d> error_ij{ residuals };

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix2d Ri = Xi.block(0, 0, 2, 2);
        Eigen::Vector2d ti{ xi(0), xi(1) };

        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix2d Rj = Xj.block(0, 0, 2, 2);
        Eigen::Vector2d tj{ xj(0), xj(1) };

        Eigen::Matrix3d Z  = PoseToTrans(measurement);
        Eigen::Matrix2d Rij = Z.block(0, 0, 2, 2);
        Eigen::Vector2d tij{ measurement(0), measurement(1) };

        Eigen::Matrix2d dRiT_dtheta;       //  derivative of Ri^T over theta
        dRiT_dtheta(0, 0) = -1 * Ri(1, 0); //  cosX -> -sinX
        dRiT_dtheta(0, 1) =  1 * Ri(0, 0); //  sinX ->  cosX
        dRiT_dtheta(1, 0) = -1 * Ri(0, 0); // -sinX -> -cosX
        dRiT_dtheta(1, 1) = -1 * Ri(1, 0); //  cosX -> -sinX

        // calcuate error & normalize error on theta
        error_ij.segment<2>(0) = Rij.transpose() * (Ri.transpose() * (tj - ti) - tij);
        error_ij(2) = xj(2) - xi(2) - measurement(2);
        if (error_ij(2) > M_PI) {
            error_ij(2) -= 2 * M_PI;
        } else if (error_ij(2) < -1 * M_PI) {
            error_ij(2) += 2 * M_PI;
        }

        error_ij = sqrt_info_matrix * error_ij;

        Eigen::Matrix3d Ai;
        Eigen::Matrix3d Bi;
        Ai.setZero();
        Ai.block(0, 0, 2, 2) = -Rij.transpose() * Ri.transpose();
        Ai.block(0, 2, 2, 1) = Rij.transpose() * dRiT_dtheta * (tj - ti);
        Ai(2, 2) = -1.0;

        Bi.setIdentity();
        Bi.block(0, 0, 2, 2) = Rij.transpose() * Ri.transpose();

        if(jacobians){
            if(jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > jacobian_xi(jacobians[0]);
                jacobian_xi = sqrt_info_matrix * Ai;
            }
            if(jacobians[1])
            {
                Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > jacobian_xj(jacobians[1]);
                jacobian_xj = sqrt_info_matrix * Bi;
            }
        }

        return true;
    }

private:
    Eigen::Vector3d measurement;
    Eigen::Matrix3d sqrt_info_matrix;
};

template <typename T>
Eigen::Matrix<T, 2, 2> RotationMatrix2D(T yaw_radians) {
  const T cos_yaw = ceres::cos(yaw_radians);
  const T sin_yaw = ceres::sin(yaw_radians);
  Eigen::Matrix<T, 2, 2> rotation;
  rotation << cos_yaw, -sin_yaw, sin_yaw, cos_yaw;
  return rotation;
}

// Normalizes the angle in radians between [-pi and pi).
template <typename T>
inline T NormalizeAngle(const T& angle_radians) {
  // Use ceres::floor because it is specialized for double and Jet types.
  T two_pi(2.0 * M_PI);
  return angle_radians -
         two_pi * ceres::floor((angle_radians + T(M_PI)) / two_pi);
}

// Computes the error term for two poses that have a relative pose measurement
// between them. Let the hat variables be the measurement.
//
// residual =  information^{1/2} * [  r_a^T * (p_b - p_a) - \hat{p_ab}   ]
//                                 [ Normalize(yaw_b - yaw_a - \hat{yaw_ab}) ]
//
// where r_a is the rotation matrix that rotates a vector represented in frame A
// into the global frame, and Normalize(*) ensures the angles are in the range
// [-pi, pi).
class PoseGraph2dErrorTerm {
 public:
  PoseGraph2dErrorTerm(double x_ab,
                       double y_ab,
                       double yaw_ab_radians,
                       const Eigen::Matrix3d& sqrt_information)
      : p_ab_(x_ab, y_ab),
        yaw_ab_radians_(yaw_ab_radians),
        sqrt_information_(sqrt_information) {}
  template <typename T>
  bool operator()(const T* const x_a,
                  const T* const y_a,
                  const T* const yaw_a,
                  const T* const x_b,
                  const T* const y_b,
                  const T* const yaw_b,
                  T* residuals_ptr) const {
    const Eigen::Matrix<T, 2, 1> p_a(*x_a, *y_a);
    const Eigen::Matrix<T, 2, 1> p_b(*x_b, *y_b);
    Eigen::Map<Eigen::Matrix<T, 3, 1>> residuals_map(residuals_ptr);
    residuals_map.template head<2>() =
        RotationMatrix2D(*yaw_a).transpose() * (p_b - p_a) - p_ab_.cast<T>();
    residuals_map(2) = NormalizeAngle(
        (*yaw_b - *yaw_a) - static_cast<T>(yaw_ab_radians_));
    // Scale the residuals by the square root information matrix to account for
    // the measurement uncertainty.
    residuals_map = sqrt_information_.template cast<T>() * residuals_map;
    return true;
  }
  static ceres::CostFunction* Create(double x_ab,
                                     double y_ab,
                                     double yaw_ab_radians,
                                     const Eigen::Matrix3d& sqrt_information) {
    return (new ceres::
                AutoDiffCostFunction<PoseGraph2dErrorTerm, 3, 1, 1, 1, 1, 1, 1>(
                    new PoseGraph2dErrorTerm(
                        x_ab, y_ab, yaw_ab_radians, sqrt_information)));
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 private:
  // The position of B relative to A in the A frame.
  const Eigen::Vector2d p_ab_;
  // The orientation of frame B relative to frame A.
  const double yaw_ab_radians_;
  // The inverse square root of the measurement covariance matrix.
  const Eigen::Matrix3d sqrt_information_;
};