#pragma once
#include "cross_section.h"

extern unsigned nMid70;
extern unsigned nHorizon;
extern unsigned nSensor;
constexpr auto SCANNING_TIME = 30U;
class Livox/* :public Cylinder, public Arch*/
{
public:
    Livox();
    Livox(float tx, float ty, float tz, float rx, float ry);
    ~Livox() { delete ray;  /*std::cout << "delete Livox\n";*/ }
    void print()const;
    unsigned get_nPts();
    static float d_to_r(float degree);
    static void visualiation(std::vector<Livox>& Sensors);

protected:
    Eigen::Vector3d* ray;

    constexpr static float scalar = 100.0f;
    float ox;
    float oy;
    float oz;
    float roll;
    float pitch;

    unsigned nPts;
    unsigned nIntersections;
   
    constexpr static unsigned scanning_time = SCANNING_TIME;
    constexpr static float rpm1 = 7294.0f;
    constexpr static float w1 = rpm1 / 60.0f * 2.0f * EIGEN_PI;
};

class Mid : public Livox
{
public:
    Mid();
    Mid(float tx, float ty, float tz, float rx, float ry, float rz, float fov);
    ~Mid() { /*std::cout << "delete Mid\n";*/ }
protected:
    float FoV;
    static unsigned nPtsPerSec;

    static float dt;
    static float rpm2;
    static float w2;
};

class Mid70 : public Mid
{
public:
    Mid70();
    Mid70(float tx, float ty, float tz, float rx, float ry);
    Mid70(float tx, float ty, float tz, float rx, float ry, float rz);
    Mid70(std::vector<float>::iterator b);
    Mid70(float* ptr);
    ~Mid70() { /*std::cout << "delete Mid70 and ray\n";*/ }
};

class Mid40 : public Mid
{
public:
    Mid40();
    Mid40(float tx, float ty, float tz, float rx, float ry);
    Mid40(std::vector<float>::iterator b);
    ~Mid40() { /*std::cout << "delete Mid40 and ray\n";*/ }
};

class Horizon : public Livox
{
public:
    Horizon();
    Horizon(float tx, float ty, float tz, float rx, float ry, float rz);
    Horizon(std::vector<float>::iterator b);
    Horizon(float* ptr);
    ~Horizon() { /*std::cout << "delete Horizon and ray\n";*/ }

private:
    float yaw;

    static float FoV_h;
    static float FoV_v;

    static unsigned nPtsPerSec;

    static float dt;
    static float rpm2;
    static float w2;
    static float rpm3;
    static float w3;
    static float n1;
    static float n2;
    static float alpha1;
    static float alpha2;
    static float delta1;
    static float delta2;
    static float rotation_ratio;
    static float translation_ratio;
};