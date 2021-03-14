#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SparseCholesky>
//#include<Eigen/IterativeLinearSolvers>

#include <iostream>
#include <ctime>

#include "utilities.hpp"

//计算整个pose-graph的误差
double ComputeError(std::vector<Eigen::Vector3d>& Vertexs,
                    std::vector<myEdge>& Edges)
{
    double sumError = 0;
    for(int i = 0; i < Edges.size();i++)
    {
        myEdge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z  = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() *  Xi.inverse() * Xj;

        Eigen::Vector3d ei = TransToPose(Ei);


        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}


/**
 * @brief CalcJacobianAndError
 *         计算jacobian矩阵和error
 * @param xi    fromIdx
 * @param xj    toIdx
 * @param z     观测值:xj相对于xi的坐标
 * @param ei    计算的误差
 * @param Ai    相对于xi的Jacobian矩阵
 * @param Bi    相对于xj的Jacobian矩阵
 */
void CalcJacobianAndError(Eigen::Vector3d& xi,Eigen::Vector3d& xj,Eigen::Vector3d& z,
                          Eigen::Vector3d& ei,Eigen::Matrix3d& Ai,Eigen::Matrix3d& Bi)
{
    //TODO--Start
    Eigen::Matrix3d Xi = PoseToTrans(xi);
    Eigen::Matrix2d Ri = Xi.block(0, 0, 2, 2);
    Eigen::Vector2d ti{ xi(0), xi(1)};

    Eigen::Matrix3d Xj = PoseToTrans(xj);
    Eigen::Matrix2d Rj = Xj.block(0, 0, 2, 2);
    Eigen::Vector2d tj{ xj(0), xj(1)};

    Eigen::Matrix3d Z  = PoseToTrans(z);
    Eigen::Matrix2d Rij = Z.block(0, 0, 2, 2);
    Eigen::Vector2d tij{ z(0), z(1)};

    Eigen::Matrix2d dRiT_dtheta;       //  derivative of Ri^T over theta
    dRiT_dtheta(0, 0) = -1 * Ri(1, 0); //  cosX -> -sinX
    dRiT_dtheta(0, 1) =  1 * Ri(0, 0); //  sinX ->  cosX
    dRiT_dtheta(1, 0) = -1 * Ri(0, 0); // -sinX -> -cosX
    dRiT_dtheta(1, 1) = -1 * Ri(1, 0); //  cosX -> -sinX

    // calcuate error & normalize error on theta
    ei.segment<2>(0) = Rij.transpose() * (Ri.transpose() * (tj - ti) - tij);
    ei(2) = xj(2) - xi(2) - z(2);
    if (ei(2) > M_PI) {
        ei(2) -= 2 * M_PI;
    } else if (ei(2) < -1 * M_PI) {
        ei(2) += 2 * M_PI;
    }

    Ai.setZero();
    Ai.block(0, 0, 2, 2) = -Rij.transpose() * Ri.transpose();
    Ai.block(0, 2, 2, 1) = Rij.transpose() * dRiT_dtheta * (tj - ti);
    Ai(2, 2) = -1.0;

    Bi.setIdentity();
    Bi.block(0, 0, 2, 2) = Rij.transpose() * Ri.transpose();

}

/**
 * @brief LinearizeAndSolve
 *        高斯牛顿方法的一次迭代．
 * @param Vertexs   图中的所有节点
 * @param Edges     图中的所有边
 * @return          位姿的增量
 */
Eigen::VectorXd  LinearizeAndSolve(std::vector<Eigen::Vector3d>& Vertexs,
                                   std::vector<myEdge>& Edges)
{
    //申请内存
    static Eigen::MatrixXd H(Vertexs.size() * 3,Vertexs.size() * 3);
    static Eigen::VectorXd b(Vertexs.size() * 3);

    H.setZero();
    b.setZero();

    //固定第一帧
    Eigen::Matrix3d I;
    I.setIdentity();
    H.block(0, 0, 3, 3) += I;

    //构造H矩阵　＆ b向量
    for(int i = 0; i < Edges.size();i++)
    {
        //提取信息
        myEdge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        //计算误差和对应的Jacobian
        Eigen::Vector3d ei;
        Eigen::Matrix3d Ai;
        Eigen::Matrix3d Bi;
        CalcJacobianAndError(xi,xj,z,ei,Ai,Bi);

        // 更新 b
        b.segment<3>(tmpEdge.xi * 3) += Ai.transpose() * infoMatrix * ei;
        b.segment<3>(tmpEdge.xj * 3) += Bi.transpose() * infoMatrix * ei;

        // 更新 H
        H.block(tmpEdge.xi * 3, tmpEdge.xi * 3, 3, 3) += Ai.transpose() * infoMatrix * Ai;
        H.block(tmpEdge.xi * 3, tmpEdge.xj * 3, 3, 3) += Ai.transpose() * infoMatrix * Bi;
        H.block(tmpEdge.xj * 3, tmpEdge.xi * 3, 3, 3) += Bi.transpose() * infoMatrix * Ai;
        H.block(tmpEdge.xj * 3, tmpEdge.xj * 3, 3, 3) += Bi.transpose() * infoMatrix * Bi;
    }

    //求解
    Eigen::VectorXd dx;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> sp_H = H.sparseView();
    solver.compute(sp_H);
    if (solver.info() != Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
    }
    dx = solver.solve(-b);
    if (solver.info() != Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
    }

    return dx;
}











