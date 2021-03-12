#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>

#include <eigen3/Eigen/Core>

//位姿-->转换矩阵
template<typename T>
Eigen::Matrix<T, 3, 3> PoseToTrans(Eigen::Matrix<T, 3, 1> x)
{
    Eigen::Matrix<T, 3, 3> trans;
    // std::cout << x(0) << std::endl;
    trans << cos(x(2)),-sin(x(2)),  x(0),
             sin(x(2)), cos(x(2)),  x(1),
                  T(0),      T(0),  T(1);

    return trans;
}


//转换矩阵－－＞位姿
template<typename T>
Eigen::Matrix<T, 3, 1> TransToPose(Eigen::Matrix<T, 3, 3> trans)
{
    Eigen::Matrix<T, 3, 1> pose;
    pose(0) = trans(0,2);
    pose(1) = trans(1,2);
    pose(2) = atan2(trans(1,0),trans(0,0));

    return pose;
}

#endif