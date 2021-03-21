#include <cmath>

#include "occupancy_mapping.h"
#include "nav_msgs/GetMap.h"
#include "sensor_msgs/PointCloud.h"
#include "sensor_msgs/PointCloud2.h"
#include "geometry_msgs/Point32.h"

/**
 * Increments all the grid cells from (x0, y0) to (x1, y1);
 * //不包含(x1,y1)
 * 2D画线算法　来进行计算两个点之间的grid cell
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 */
std::vector<GridIndex> TraceLine(int x0, int y0, int x1, int y1)
{
    GridIndex tmpIndex;
    std::vector<GridIndex> gridIndexVector;

    bool steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep)
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    if (x0 > x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int deltaX = x1 - x0;
    int deltaY = abs(y1 - y0);
    int error = 0;
    int ystep;
    int y = y0;

    if (y0 < y1)
    {
        ystep = 1;
    }
    else
    {
        ystep = -1;
    }

    int pointX;
    int pointY;
    for (int x = x0; x <= x1; x++)
    {
        if (steep)
        {
            pointX = y;
            pointY = x;
        }
        else
        {
            pointX = x;
            pointY = y;
        }

        error += deltaY;

        if (2 * error >= deltaX)
        {
            y += ystep;
            error -= deltaX;
        }

        //不包含最后一个点．
        if (pointX == x1 && pointY == y1)
            continue;

        //保存所有的点
        tmpIndex.SetIndex(pointX, pointY);

        gridIndexVector.push_back(tmpIndex);
    }

    return gridIndexVector;
}

void SetMapParams(void)
{
    mapParams.width = 1000;
    mapParams.height = 1000;
    mapParams.resolution = 0.05;

    //每次被击中的log变化值，覆盖栅格建图算法需要的参数
    mapParams.log_free = -1;
    mapParams.log_occ = 3;

    //每个栅格的最大最小值．
    mapParams.log_max = 100.0;
    mapParams.log_min = 0.0;

    mapParams.origin_x = 0.0;
    mapParams.origin_y = 0.0;

    //地图的原点，在地图的正中间
    mapParams.offset_x = 500;
    mapParams.offset_y = 500;

    pMap = new unsigned char[mapParams.width * mapParams.height];

    //计数建图算法需要的参数
    //每个栅格被激光击中的次数
    pMapHits = new unsigned long[mapParams.width * mapParams.height];
    //每个栅格被激光通过的次数
    pMapMisses = new unsigned long[mapParams.width * mapParams.height];

    //TSDF建图算法需要的参数
    pMapW = new unsigned long[mapParams.width * mapParams.height];
    pMapTSDF = new double[mapParams.width * mapParams.height];

    //初始化
    for (int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        pMap[i] = 50;
        pMapHits[i] = 0;
        pMapMisses[i] = 0;
        pMapW[i] = 0;
        pMapTSDF[i] = -1;
    }
}

//从世界坐标系转换到栅格坐标系
GridIndex ConvertWorld2GridIndex(double x, double y)
{
    GridIndex index;

    index.x = std::ceil((x - mapParams.origin_x) / mapParams.resolution) + mapParams.offset_x;
    index.y = std::ceil((y - mapParams.origin_y) / mapParams.resolution) + mapParams.offset_y;

    return index;
}

int GridIndexToLinearIndex(const GridIndex& index)
{
    int linear_index;
    linear_index = index.y + index.x * mapParams.width;
}

//判断index是否有效
bool isValidGridIndex(const GridIndex& index)
{
    if (index.x >= 0 && index.x < mapParams.width && index.y >= 0 && index.y < mapParams.height)
        return true;

    return false;
}

void DestoryMap()
{
    if (pMap != NULL)
        delete pMap;
}

void OccupancyGridUpdate(const GridIndex& end_grid_index, const GridIndex& robot_index) 
{
    auto updateGrid = [](const GridIndex& index, double update_value) -> void {
        // convert pMap to log-odd then update
        int linear_grid_index = GridIndexToLinearIndex(index);
        double px = (double) pMap[linear_grid_index] / 100.0;
        double lx = std::log(px / (1 - px));
        lx += (update_value);
        
        // convert log-odd back to pMap
        pMap[linear_grid_index] = 100 * (1.0 - (1.0 / (1 + std::exp(lx))));
        pMap[linear_grid_index] = std::max(mapParams.log_min, std::min(mapParams.log_max, (double)pMap[linear_grid_index]));
    };

    // update end grid occupied by lidar beam with log-odds of occupation
    if (isValidGridIndex(end_grid_index)) {
        int linear_grid_index = GridIndexToLinearIndex(end_grid_index);
        updateGrid(end_grid_index, mapParams.log_occ);
    }

    // update grid that are passed through by lidar beam with log-odds of free
    std::vector<GridIndex> free_grid_indices = TraceLine(robot_index.x, robot_index.y, end_grid_index.x, end_grid_index.y);
    for (const GridIndex& grid_index: free_grid_indices) {
        if (isValidGridIndex(grid_index)) {
            updateGrid(grid_index, mapParams.log_free);
        }
    }
}

void CountMappingUpdate(const GridIndex& end_grid_index, const GridIndex& robot_index)
{
    if (isValidGridIndex(end_grid_index)) {
        int linear_grid_index = GridIndexToLinearIndex(end_grid_index);
        pMapHits[linear_grid_index]++;
    }

    // update grid that are passed through by lidar beam with log-odds of free
    std::vector<GridIndex> free_grid_indices = TraceLine(robot_index.x, robot_index.y, end_grid_index.x, end_grid_index.y);
    for (const GridIndex& grid_index: free_grid_indices) {
        if (isValidGridIndex(grid_index)) {
            int linear_grid_index = GridIndexToLinearIndex(grid_index);
            pMapMisses[linear_grid_index]++;
        }
    }
}

void TSDFMappingUpdate(const GridIndex& end_grid_index, const GridIndex& robot_index, double laser_dist)
{  
    // number of grids surronding by the end grid index that we update
    static const int nums_surronding_grids = 10;
    double robot_laser_x = end_grid_index.x - robot_index.x;
    double robot_laser_y = end_grid_index.y - robot_index.y;
    
    for (int i = -nums_surronding_grids; i <= nums_surronding_grids; ++i) {
        GridIndex index = end_grid_index;
        double dx = 0;
        double dy = 0;
        
        if (robot_laser_x == 0) {
            // dx = 0;
            dy = i;
        } else if (robot_laser_y == 0) {
            // dy = i;
            dx = i;
        } else {
            dx = i;
            dy = i * std::fabs(robot_laser_y / robot_laser_x);
        }

        index.x += dx;
        index.y += dy;

        if (!isValidGridIndex(index)) {
            continue;
        }

        double grid_robot_dist = std::sqrt(std::pow(index.x - robot_index.x, 2) + std::pow(index.y - robot_index.y, 2));
        grid_robot_dist *= mapParams.resolution;
        double sdf = laser_dist - grid_robot_dist;
        double tsdf = std::min(1.0, std::max(-1.0, sdf / 0.1));

        int linear_grid_index = GridIndexToLinearIndex(index);

        // use 1.0 as weight
        pMapTSDF[linear_grid_index] = (pMapW[linear_grid_index] * pMapTSDF[linear_grid_index] + 1.0 * tsdf) / (1.0 + pMapW[linear_grid_index]);
        pMapW[linear_grid_index] += 1.0;
    }
}

//
void OccupancyMapping(std::vector<GeneralLaserScan> &scans, std::vector<Eigen::Vector3d> &robot_poses)
{
    std::cout << "开始建图，请稍后..." << std::endl;
    //枚举所有的激光雷达数据
    for (int i = 0; i < scans.size(); i++)
    {
        GeneralLaserScan scan = scans[i];
        Eigen::Vector3d robotPose = robot_poses[i];

        //机器人的下标
        GridIndex robotIndex = ConvertWorld2GridIndex(robotPose(0), robotPose(1));

        for (int id = 0; id < scan.range_readings.size(); id++)
        {
            double dist = scan.range_readings[id];
            double angle = -scan.angle_readings[id]; // 激光雷达逆时针转，角度取反

            if (std::isinf(dist) || std::isnan(dist))
                continue;

            //计算得到该激光点的世界坐标系的坐标
            double theta = -robotPose(2); // 激光雷达逆时针转，角度取反
            double laser_x = dist * cos(angle);
            double laser_y = dist * sin(angle);

            double world_x = cos(theta) * laser_x - sin(theta) * laser_y + robotPose(0);
            double world_y = sin(theta) * laser_x + cos(theta) * laser_y + robotPose(1);

            //start of TODO 对对应的map的cell信息进行更新．（1,2,3题内容）
            GridIndex end_grid_index = ConvertWorld2GridIndex(world_x, world_y);
            // OccupancyGridUpdate(end_grid_index, robotIndex);
            // CountMappingUpdate(end_grid_index, robotIndex);
            TSDFMappingUpdate(end_grid_index, robotIndex, dist);
            //end of TODO
        }
    }
    //start of TODO 通过计数建图算法或TSDF算法对栅格进行更新（2,3题内容）
    // use CountMapping to update occupation probability
    // for (size_t i = 0; i < mapParams.height * mapParams.width; ++i) {
    //     // Count mapping final update
    //     if (pMapMisses[i] + pMapHits[i] != 0) {
    //         pMap[i] = 100 * static_cast<int>(static_cast<double>(pMapHits[i]) /static_cast<double>(pMapMisses[i] + pMapHits[i]));
    //     }
    // }

    // use pMap to display fielf of TSDF
    for (size_t i = 0; i < mapParams.height * mapParams.width; ++i) {
        // Count mapping final update
        if (pMapTSDF[i] == -1) {
            continue;
        }
        if (pMapTSDF[i] == 1) {
            pMap[i] = 0.0;
            std::cout << "Found free grid.\n";
            continue;
        }
        pMap[i] = 100 * ((pMapTSDF[i] + 1.0) / 2.0);
        // pMap[i] = pMapTSDF[i] < 0 ? 80 : 20;
    }

    // use TSDF to update occupation probability
    for (size_t i = 0; i < mapParams.height; ++i) {
        double prev_tsdf = pMapTSDF[i * mapParams.width];
        for (size_t j = 1; j < mapParams.width; ++j) {
            double curr_tsdf = pMapTSDF[i * mapParams.width + j];

            // mark pMap as free/unknown/occupied based on sign of TSDF value:
            // 1. negative -> unknown
            // 2. positive -> free
            // 3. exact zero(rarely happens!) -> occupied
            if (curr_tsdf > 0) {
                pMap[i * mapParams.width + j] = 0;
            } else if (curr_tsdf < 0) {
                pMapW[i * mapParams.width + j] = -1;
            } else {
                pMapW[i * mapParams.width + j] = 100;
            }

            // linear interpolation to find the grid where 0 tsdf lies on
            if (prev_tsdf * curr_tsdf < 0) {
                if (std::fabs(prev_tsdf) > std::fabs(curr_tsdf)) {
                    pMap[i * mapParams.width + j] = 100;
                } else {
                    pMap[i * mapParams.width + j - 1] = 100;
                }
            }

            prev_tsdf = curr_tsdf;
        }
    }
    //end of TODO
    std::cout << "建图完毕" << std::endl;
}

//发布地图．
void PublishMap(ros::Publisher &map_pub)
{
    nav_msgs::OccupancyGrid rosMap;

    rosMap.info.resolution = mapParams.resolution;
    rosMap.info.origin.position.x = 0.0;
    rosMap.info.origin.position.y = 0.0;
    rosMap.info.origin.position.z = 0.0;
    rosMap.info.origin.orientation.x = 0.0;
    rosMap.info.origin.orientation.y = 0.0;
    rosMap.info.origin.orientation.z = 0.0;
    rosMap.info.origin.orientation.w = 1.0;

    rosMap.info.origin.position.x = mapParams.origin_x;
    rosMap.info.origin.position.y = mapParams.origin_y;
    rosMap.info.width = mapParams.width;
    rosMap.info.height = mapParams.height;
    rosMap.data.resize(rosMap.info.width * rosMap.info.height);

    //0~100
    int cnt0, cnt1, cnt2;
    cnt0 = cnt1 = cnt2 = 100;
    for (int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        if (pMap[i] == 50)
        {
            rosMap.data[i] = -1.0;
        }
        else
        {

            rosMap.data[i] = pMap[i];
        }
    }

    rosMap.header.stamp = ros::Time::now();
    rosMap.header.frame_id = "map";

    map_pub.publish(rosMap);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "OccupancyMapping");

    ros::NodeHandle nodeHandler;

    ros::Publisher mapPub = nodeHandler.advertise<nav_msgs::OccupancyGrid>("laser_map", 1, true);

    std::vector<Eigen::Vector3d> robotPoses;
    std::vector<GeneralLaserScan> generalLaserScans;

    std::string basePath = "/home/xt/Projects/LidarSLAM-learning-record/hw7/OccupancyMappingProject/src/data";

    std::string posePath = basePath + "/pose.txt";
    std::string anglePath = basePath + "/scanAngles.txt";
    std::string scanPath = basePath + "/ranges.txt";

    //读取数据
    ReadPoseInformation(posePath, robotPoses);

    ReadLaserScanInformation(anglePath,
                             scanPath,
                             generalLaserScans);

    //设置地图信息
    SetMapParams();

    OccupancyMapping(generalLaserScans, robotPoses);

    PublishMap(mapPub);

    ros::spin();

    DestoryMap();

    std::cout << "Release Memory!!" << std::endl;
}
