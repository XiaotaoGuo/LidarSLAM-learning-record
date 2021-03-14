#ifndef G2O_METHOD_H
#define G2O_METHOD_H

#include <g2o/core/g2o_core_api.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/sparse_optimizer_terminate_action.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>


#include "ls_slam/my_types.h"

namespace g2o_method
{

void solveProblems(std::vector<myVertex>& vertexes, std::vector<myEdge>& edges);

class MyG2OVertex : public g2o::BaseVertex<3, Eigen::Vector3d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual void setToOriginImpl() override {
        _estimate << 0, 0, 0;
    }

    virtual void oplusImpl(const double *update) override {
        _estimate += Eigen::Vector3d(update);
    }

    virtual bool read(std::istream &in) {}
    virtual bool write(std::ostream &out) const {}

};

class MyG2OEdge : public g2o::BaseBinaryEdge<3, Eigen::Vector3d, MyG2OVertex, MyG2OVertex> 
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    
    virtual void computeError() override;
    virtual void linearizeOplus() override;  

    virtual bool read(std::istream& is) {}
    virtual bool write(std::ostream& os) const {}  
};
    
} // namespace g2o_method


#endif