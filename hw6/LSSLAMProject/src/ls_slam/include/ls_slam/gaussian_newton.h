#ifndef GAUSSIAN_NEWTON_H
#define GAUSSIAN_NEWTON_H

#include <vector>
#include <eigen3/Eigen/Core>

#include "my_types.h"

Eigen::VectorXd  LinearizeAndSolve(std::vector<myVertex>& Vertexs,
                                   std::vector<myEdge>& Edges);

double ComputeError(std::vector<myVertex>& Vertexs,
                    std::vector<myEdge>& Edges);








#endif
