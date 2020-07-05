#include <map.h>
#include "gaussian_newton_method.h"

const double GN_PI = 3.1415926;

//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

Eigen::Matrix3d GN_V2T(Eigen::Vector3d vec)
{
    Eigen::Matrix3d T;
    T  << cos(vec(2)),-sin(vec(2)),vec(0),
            sin(vec(2)), cos(vec(2)),vec(1),
            0,           0,     1;

    return T;
}

//对某一个点进行转换．
Eigen::Vector2d GN_TransPoint(Eigen::Vector2d pt,Eigen::Matrix3d T)
{
    Eigen::Vector3d tmp_pt(pt(0),pt(1),1);
    tmp_pt = T * tmp_pt;
    return Eigen::Vector2d(tmp_pt(0),tmp_pt(1));
}



//用激光雷达数据创建势场．
map_t* CreateMapFromLaserPoints(Eigen::Vector3d map_origin_pt,
                                std::vector<Eigen::Vector2d> laser_pts,
                                double resolution)
{
    map_t* map = map_alloc();

    map->origin_x = map_origin_pt(0);
    map->origin_y = map_origin_pt(1);
    map->resolution = resolution;

    //固定大小的地图，必要时可以扩大．
    map->size_x = 10000;
    map->size_y = 10000;

    map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);

    //高斯平滑的sigma－－固定死
    map->likelihood_sigma = 0.5;

    Eigen::Matrix3d Trans = GN_V2T(map_origin_pt);

    //设置障碍物
    for(int i = 0; i < laser_pts.size();i++)
    {
        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);

        int cell_x,cell_y;
        cell_x = MAP_GXWX(map,tmp_pt(0));
        cell_y = MAP_GYWY(map,tmp_pt(1));

        map->cells[MAP_INDEX(map,cell_x,cell_y)].occ_state = CELL_STATUS_OCC;
    }

    //进行障碍物的膨胀--最大距离固定死．
    map_update_cspace(map,0.5);

    return map;
}


/**
 * @brief InterpMapValueWithDerivatives
 * 在地图上的进行插值，得到coords处的势场值和对应的关于位置的梯度．
 * 返回值为Eigen::Vector3d ans
 * ans(0)表示势场值
 * ans(1:2)表示梯度
 * @param map
 * @param coords
 * @return
 */
Eigen::Vector3d InterpMapValueWithDerivatives(map_t* map,Eigen::Vector2d& coords)
{
    Eigen::Vector3d ans;

    /// calculate coordinates on the given map (noticed: might not be exact interger)
    /// remember the origin of map locates in center of the map thus need to plus half of size to get the real index
    double index_x = (coords(0) - map->origin_x) / map->resolution + map->size_x / 2;
    double index_y = (coords(1) - map->origin_y) / map->resolution + map->size_y / 2;
    int16_t index_x0 = floor(index_x);
    int16_t index_y0 = floor(index_y);

    /// calculate u, v using four nearest points
    double u, v;
    u = index_x - index_x0;
    v = index_y - index_y0;

    /// calcualte scores for four nearest points
    double z1 = map->cells[MAP_INDEX(map, index_x0, index_y0)].score;
    double z2 = map->cells[MAP_INDEX(map, index_x0 + 1, index_y0)].score;
    double z3 = map->cells[MAP_INDEX(map, index_x0 + 1, index_y0 + 1)].score;
    double z4 = map->cells[MAP_INDEX(map, index_x0, index_y0 + 1)].score;

    /// score of given coordinate in the map
    ans(0) = (1 - u) * (1 - v) * z1 + u * (1 - v) * z2 + u * v * z3 + (1 - u) * v * z4;

    /// gradient
    /// noticed: need to remove scale influence by using resolution
    ans(1) = (v * (z3 - z4) + (1 - v) * (z2 - z1)) / map->resolution;
    ans(2) = (u * (z3 - z2) + (1 - u) * (z4 - z1)) / map->resolution;

    return ans;
}


/**
 * @brief ComputeCompleteHessianAndb
 * 计算H*dx = b中的H和b
 * @param map
 * @param now_pose
 * @param laser_pts
 * @param H
 * @param b
 */
void ComputeHessianAndb(map_t* map, Eigen::Vector3d now_pose,
                        std::vector<Eigen::Vector2d>& laser_pts,
                        Eigen::Matrix3d& H, Eigen::Vector3d& b)
{
    H = Eigen::Matrix3d::Zero();
    b = Eigen::Vector3d::Zero();

    Eigen::Matrix3d T = GN_V2T(now_pose);

    for (Eigen::Vector2d pt: laser_pts) {

        /// compute laser pt pose globally
        Eigen::Vector2d pt_pose = GN_TransPoint(pt, T);

        Eigen::Matrix<double, 2, 3> ds;
        ds << 1, 0, -sin(now_pose(2)) * pt(0) - cos(now_pose(2)) * pt(1),
                0, 1, cos(now_pose(2)) * pt(0) - sin(now_pose(2)) * pt(1);

        /// compute score & grident of that point in map
        Eigen::Vector3d score_gradient = InterpMapValueWithDerivatives(map, pt_pose);
        Eigen::Vector2d gradient(score_gradient(1), score_gradient(2));
        double score = score_gradient(0);

        /// noticed the dimension of J should 1 x 3
        Eigen::RowVector3d J = gradient.transpose() * ds;
        H += J.transpose() * J;
        b += J.transpose() * (1 - score);
    }

}


/**
 * @brief GaussianNewtonOptimization
 * 进行高斯牛顿优化．
 * @param map
 * @param init_pose
 * @param laser_pts
 */
void GaussianNewtonOptimization(map_t*map,Eigen::Vector3d& init_pose,std::vector<Eigen::Vector2d>& laser_pts)
{
    int maxIteration = 20;
    Eigen::Vector3d now_pose = init_pose;
    Eigen::Matrix3d H;
    Eigen::Vector3d b;

    for(int i = 0; i < maxIteration;i++)
    {

        ComputeHessianAndb(map, now_pose, laser_pts, H, b);

        Eigen::Vector3d delta_x = H.colPivHouseholderQr().solve(b);

        now_pose += delta_x;
    }
    init_pose = now_pose;

}
