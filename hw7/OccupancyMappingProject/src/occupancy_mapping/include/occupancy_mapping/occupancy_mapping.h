#ifndef OCCUPANCY_MAPPING_H
#define OCCUPANCY_MAPPING_H

#include <iostream>
#include <vector>

#include <ros/ros.h>

#include <eigen3/Eigen/Core>

#include "readfile.h"

typedef struct gridindex_
{
    int x;
    int y;

    void SetIndex(int x_, int y_)
    {
        x = x_;
        y = y_;
    }
} GridIndex;

typedef struct map_params
{
    double log_occ, log_free;
    double log_max, log_min;
    double resolution;
    double origin_x, origin_y;
    int height, width;
    int offset_x, offset_y;
} MapParams;

MapParams mapParams;

// 栅格占用概率，初始值为 50
unsigned char *pMap;

// 栅格被击中和穿过的次数
unsigned long *pMapHits;
unsigned long *pMapMisses;

// 栅格的权重以及 TSDF 函数值
unsigned long *pMapW;
double *pMapTSDF;

#endif
