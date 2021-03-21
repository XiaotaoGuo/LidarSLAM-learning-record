#include "ls_slam/auto_solvers/gtsam/gtsam_method.h"

void gtsam_method::solveProblems(std::vector<myVertex>& vertexes, std::vector<myEdge>& edges)
{
    gtsam::NonlinearFactorGraph graph;

    gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector3(0.2, 0.2, 0.1));
    graph.addPrior(0, gtsam::Pose2(0, 0, 0), priorNoise);

    for (const auto& edge: edges) {
        gtsam::Matrix sqrtInfo = edge.infoMatrix.array().sqrt();
        gtsam::Pose2 measurement(edge.measurement(0), edge.measurement(1), edge.measurement(2));
        graph.emplace_shared<gtsam::BetweenFactor<gtsam::Pose2>>(edge.xi, edge.xj,\
                                                                measurement, \
                                                                gtsam::noiseModel::Gaussian::SqrtInformation(sqrtInfo));
    }

    gtsam::Values initialEstimate;
    for (size_t i = 0; i < vertexes.size(); i++) {
        gtsam::Pose2 v(vertexes[i](0), vertexes[i](1), vertexes[i](2));
        initialEstimate.insert(i, v);
    }

    gtsam::GaussNewtonParams parameters;
    parameters.relativeErrorTol = 1e-6;
    parameters.maxIterations = 10;
    parameters.setVerbosity("SILENT");
    gtsam::GaussNewtonOptimizer optimizer(graph, initialEstimate, parameters);
    gtsam::Values result = optimizer.optimize();

    for (const auto& value: result) {
        gtsam::Pose2 vp = value.value.cast<gtsam::Pose2>();
        vertexes[value.key](0) = vp.x();
        vertexes[value.key](1) = vp.y();
        vertexes[value.key](2) = vp.theta();
    }
}