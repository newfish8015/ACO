#pragma once
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>
#include <string>
#include <random>
#include <chrono>

#include <ctime>
#include <cmath>
#include <cstdlib>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pcl/common/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>

#include <pcl/visualization/pcl_visualizer.h>
#include <omp.h>

extern unsigned THREADS;
extern char flag_cross_section;
extern pcl::PointCloud<pcl::PointXYZ> cloud;
extern pcl::PointCloud<pcl::PointXYZ> cloud_cube;
extern float t, t_s, t_p;

//橫切面
enum cross_section
{
    choose_cylinder = 1,
    choose_arch = 2,
    choose_horseshoe = 3,
    choose_square = 4
};

class Cross_Section
{
public:
    Cross_Section();
    ~Cross_Section() { delete cube; }
    static std::vector<unsigned>::iterator check_continue(std::vector<unsigned>::iterator i, std::vector<unsigned>&);

protected:
    constexpr static float eps = 1e-5f;
    constexpr static float seg = 50.0;
    constexpr static float start = 0.0;
    constexpr static float end = start + seg;
    constexpr static float precision = 0.1f;
    constexpr static float scale = 1.0f / precision;//一公尺有幾個方塊

    static unsigned cube_x;
    static unsigned cube_y;
    static unsigned cube_z;
    static unsigned nCube;
    static unsigned* cube;
};

class Cylinder : public Cross_Section
{
public:
    Cylinder();
    ~Cylinder() {}
    //static void cylinder(Eigen::Vector3d*& v, unsigned length);
    static unsigned cylinder(std::vector<Eigen::Vector3d>& v, unsigned length);
    static float statistics(bool print = false);
    constexpr static float r = 6.0f;
protected:
    constexpr static float h = 0.0f;
};

class Arch : public Cross_Section
{
public:
    Arch();
    ~Arch() {}
    //static void arch(Eigen::Vector3d*& v, unsigned length);
    static unsigned arch(std::vector<Eigen::Vector3d>& v, unsigned length);
    static float statistics(bool print = false);
    constexpr static float l = 6.0f;
protected:
    constexpr static float alpha = 45.0f / 180.0f * EIGEN_PI;
    constexpr static float r = 3.0f / 5.0f * l;
    static float a;
    static float b;
    static float R;
    static float ix;
    static float iy;
    static float f;
};

class Horseshoe : public Cross_Section
{
public:
    Horseshoe();
    ~Horseshoe() {}
    //static void arch(Eigen::Vector3d*& v, unsigned length);
    static unsigned horseshoe(std::vector<Eigen::Vector3d>& v, unsigned length);
    static float statistics(bool print = false);
    constexpr static float l = 9.0f;
    constexpr static float r = 5.0f;
protected:
   
    constexpr static float h = l - r;
};

class Square : public Cross_Section
{
public:
    Square();
    ~Square() {}
    //static void arch(Eigen::Vector3d*& v, unsigned length);
    static unsigned square(std::vector<Eigen::Vector3d>& v, unsigned length);
    static float statistics(bool print = false);
    constexpr static float l = 5.0f;
    constexpr static float r = 2.5f;
protected:

};