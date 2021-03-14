cmake_minimum_required(VERSION 3.1.0)
project(ls_slam)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
)

find_package(Eigen3 REQUIRED)

message("USE_CERES=${USE_CERES}")
if( USE_CERES )
  add_definitions(-DUSE_CERES)  
  find_package(Ceres REQUIRED)
endif()

message("USE_CERES=${USE_G2O}")
if (USE_G2O)
  add_definitions(-DUSE_G2O)
  find_package(G2O REQUIRED)
  message(STATUS "Find G2O: ${G2O_INCLUDE_DIR}")
endif()

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES ls_slam
#  CATKIN_DEPENDS roscpp rospy std_msgs
#  DEPENDS system_lib
)

###########
## Build ##
###########
## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  include/ls_slam
  ${catkin_INCLUDE_DIRS}  ${EIGEN_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${G2O_INCLUDE_DIR}
)

add_executable(ls_slam src/main.cpp src/readfile.cpp src/gaussian_newton.cpp)
target_link_libraries(ls_slam ${catkin_LIBRARIES})

if ( USE_CERES )
  target_link_libraries(ls_slam ceres)
  target_link_libraries(ls_slam ${CSPARSE_LIBRARY} ceres)
endif()

if ( USE_G2O )
  target_sources(ls_slam PRIVATE src/g2o/g2o_method.cpp)
  include_directories(${G2O_INCLUDE_DIR})
  target_link_libraries(ls_slam ${G2O_CORE_LIBRARY} ${G2O_TYPES_SLAM2D} ${G2O_SOLVER_EIGEN} ${G2O_STUFF_LIBRARY})
endif()



#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_ls_slam.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
