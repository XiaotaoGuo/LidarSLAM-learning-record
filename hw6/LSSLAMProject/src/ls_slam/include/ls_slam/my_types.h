

#ifndef MY_TYPES_HPP
#define MY_TYPES_HPP

#include <eigen3/Eigen/Core>

/**
 * Custom Edge and Vertex 
 */
struct myEdge
{
  int xi,xj;
  Eigen::Vector3d measurement;
  Eigen::Matrix3d infoMatrix;
};

typedef Eigen::Vector3d myVertex;

#endif