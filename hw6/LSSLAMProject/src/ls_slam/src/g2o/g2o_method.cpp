#include "ls_slam/auto_solvers/g2o/g2o_method.h"
#include "ls_slam/utilities.hpp"

void g2o_method::solveProblems(std::vector<myVertex>& vertexes, std::vector<myEdge>& edges)
{
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<3,3>> BlockSolverType;
    typedef g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType> LinearSolverType;

    g2o::SparseOptimizer optimizer;
    auto linearSolver = g2o::make_unique<LinearSolverType>();
    linearSolver->setBlockOrdering(false);
    g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(
        g2o::make_unique<BlockSolverType>(std::move(linearSolver)));

    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(true);

    for (size_t i = 0; i < vertexes.size(); ++i) {
        const myVertex& v_vec = vertexes[i];
        g2o::VertexSE2* v = new g2o::VertexSE2;
        g2o::SE2 se2;
        se2.fromVector(v_vec);
        v->setId(i);
        v->setEstimate(se2);
        optimizer.addVertex(v);

        // const myVertex& v_vec = vertexes[i];
        // MyG2OVertex* v = new MyG2OVertex();
        // v->setId(i);
        // v->setEstimate(v_vec);
        // optimizer.addVertex(v);
    }

    for (const auto& edge: edges) {
        g2o::EdgeSE2* e = new g2o::EdgeSE2;
        e->vertices()[0] = optimizer.vertex(edge.xi);
        e->vertices()[1] = optimizer.vertex(edge.xj);
        g2o::SE2 se2;
        se2.fromVector(edge.measurement);
        e->setMeasurement(se2);
        e->setInformation(edge.infoMatrix);
        optimizer.addEdge(e);

        // MyG2OEdge* e = new MyG2OEdge;
        // e->vertices()[0] = optimizer.vertex(edge.xi);
        // e->vertices()[1] = optimizer.vertex(edge.xj);
        // e->setMeasurement(edge.measurement);
        // e->setInformation(edge.infoMatrix);
        // optimizer.addEdge(e);

    }

    g2o::VertexSE2* firstVertex = dynamic_cast<g2o::VertexSE2*>(optimizer.vertex(0));
    // MyG2OVertex* firstVertex = dynamic_cast<MyG2OVertex*>(optimizer.vertex(0));
    firstVertex->setFixed(true);

    
    g2o::SparseOptimizerTerminateAction* terminateAction = new g2o::SparseOptimizerTerminateAction;
    terminateAction->setGainThreshold(1e-6);
    terminateAction->setMaxIterations(10);
    optimizer.addPostIterationAction(terminateAction);
    
    optimizer.initializeOptimization();
    optimizer.optimize(10);

    for (size_t i = 0; i < vertexes.size(); ++i) {
        dynamic_cast<g2o::VertexSE2*>(optimizer.vertex(i))->getEstimateData(vertexes[i].data());
        // dynamic_cast<MyG2OVertex*>(optimizer.vertex(i))->getEstimateData(vertexes[i].data());
        // vertexes[i] = dynamic_cast<MyG2OVertex*>(optimizer.vertex(i))->estimate();
    }
}

void g2o_method::MyG2OEdge::computeError() {

        MyG2OVertex* vtx1 = static_cast<MyG2OVertex*>(_vertices[0]);
        MyG2OVertex* vtx2 = static_cast<MyG2OVertex*>(_vertices[1]);
        Eigen::Vector3d measurement{ _measurement };
        Eigen::Vector3d v1 = vtx1->estimate();
        Eigen::Vector3d v2 = vtx2->estimate();

        // use sqrt(info) * Z^(-1) * X1^(-1) & X2, involving two inverse operation
        // Eigen::Matrix<T, 3, 3> X1 = PoseToTrans(v1_mat);
        // Eigen::Matrix<T, 3, 3> X2 = PoseToTrans(v2_mat);
        // Eigen::Matrix<T, 3, 3> Z = PoseToTrans(m);

        // Eigen::Matrix<T, 3, 3> E = sqrt_info_matrix.template cast<T>() * Z.inverse() * (X1.inverse() * X2);
        // error = TransToPose(E);

        // calculate error from translation and rotation respectively and combine them together
        Eigen::Matrix3d X1 = PoseToTrans(v1);
        Eigen::Matrix3d X2 = PoseToTrans(v2);
        Eigen::Matrix3d Z = PoseToTrans(measurement);
        
        Eigen::Matrix2d Ri = X1.block(0, 0, 2, 2);
        Eigen::Matrix2d Rj = X2.block(0, 0, 2, 2);
        Eigen::Matrix2d Rij = Z.block(0, 0, 2, 2);

        Eigen::Vector2d ti{ v1(0), v1(1) };
        Eigen::Vector2d tj{ v2(0), v2(1) };
        Eigen::Vector2d tij{ measurement(0), measurement(1) };

        Eigen::Matrix2d dRiT_dtheta;       //  derivative of Ri^T over theta
        dRiT_dtheta(0, 0) = -1 * Ri(1, 0); //  cosX -> -sinX
        dRiT_dtheta(0, 1) =  1 * Ri(0, 0); //  sinX ->  cosX
        dRiT_dtheta(1, 0) = -1 * Ri(0, 0); // -sinX -> -cosX
        dRiT_dtheta(1, 1) = -1 * Ri(1, 0); //  cosX -> -sinX

        // calcuate error & normalize error on theta
        _error.segment<2>(0) = Rij.transpose() * (Ri.transpose() * (tj - ti) - tij);
        _error(2) = v2(2) - v1(2) - measurement(2);
        if (_error(2) > M_PI) {
            _error(2) -= 2 * M_PI;
        } else if (_error(2) < -1 * M_PI) {
            _error(2) += 2 * M_PI;
        }
}

void g2o_method::MyG2OEdge::linearizeOplus() {

    MyG2OVertex* vtx1 = static_cast<MyG2OVertex*>(_vertices[0]);
        MyG2OVertex* vtx2 = static_cast<MyG2OVertex*>(_vertices[1]);
        Eigen::Vector3d measurement{ _measurement };
        Eigen::Vector3d v1 = vtx1->estimate();
        Eigen::Vector3d v2 = vtx2->estimate();

        // use sqrt(info) * Z^(-1) * X1^(-1) & X2, involving two inverse operation
        // Eigen::Matrix<T, 3, 3> X1 = PoseToTrans(v1_mat);
        // Eigen::Matrix<T, 3, 3> X2 = PoseToTrans(v2_mat);
        // Eigen::Matrix<T, 3, 3> Z = PoseToTrans(m);

        // Eigen::Matrix<T, 3, 3> E = sqrt_info_matrix.template cast<T>() * Z.inverse() * (X1.inverse() * X2);
        // error = TransToPose(E);

        // calculate error from translation and rotation respectively and combine them together
        Eigen::Matrix3d X1 = PoseToTrans(v1);
        Eigen::Matrix3d X2 = PoseToTrans(v2);
        Eigen::Matrix3d Z = PoseToTrans(measurement);
        
        Eigen::Matrix2d Ri = X1.block(0, 0, 2, 2);
        Eigen::Matrix2d Rj = X2.block(0, 0, 2, 2);
        Eigen::Matrix2d Rij = Z.block(0, 0, 2, 2);

        Eigen::Vector2d ti{ v1(0), v1(1) };
        Eigen::Vector2d tj{ v2(0), v2(1) };
        Eigen::Vector2d tij{ measurement(0), measurement(1) };

        Eigen::Matrix2d dRiT_dtheta;       //  derivative of Ri^T over theta
        dRiT_dtheta(0, 0) = -1 * Ri(1, 0); //  cosX -> -sinX
        dRiT_dtheta(0, 1) =  1 * Ri(0, 0); //  sinX ->  cosX
        dRiT_dtheta(1, 0) = -1 * Ri(0, 0); // -sinX -> -cosX
        dRiT_dtheta(1, 1) = -1 * Ri(1, 0); //  cosX -> -sinX

        Eigen::Matrix3d Ai;
        Eigen::Matrix3d Bi;
        Ai.setZero();
        Ai.block(0, 0, 2, 2) = -Rij.transpose() * Ri.transpose();
        Ai.block(0, 2, 2, 1) = Rij.transpose() * dRiT_dtheta * (tj - ti);
        Ai(2, 2) = -1.0;

        Bi.setIdentity();
        Bi.block(0, 0, 2, 2) = Rij.transpose() * Ri.transpose();

        _jacobianOplusXi = Ai;
        _jacobianOplusXj = Bi;
}