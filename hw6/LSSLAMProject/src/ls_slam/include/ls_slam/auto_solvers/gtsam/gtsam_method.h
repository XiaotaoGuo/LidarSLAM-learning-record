#include <gtsam/geometry/Pose2.h>
#include <gtsam/inference/Key.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/Values.h>

#include "ls_slam/my_types.h"

namespace gtsam_method
{
    void solveProblems(std::vector<myVertex>& vertexes, std::vector<myEdge>& edges);   
} // namespace gtsam_method