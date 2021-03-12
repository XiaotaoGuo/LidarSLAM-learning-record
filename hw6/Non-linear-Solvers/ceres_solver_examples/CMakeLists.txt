cmake_minimum_required(VERSION 3.5)

project(helloworld)

find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)

message(STATUS ${EIGEN3_INCLUDE_DIR})
include_directories( ${EIGEN3_INCLUDE_DIR} )

# helloworld
add_executable(helloworld ceres_hello_world.cpp)
target_link_libraries(helloworld ceres)