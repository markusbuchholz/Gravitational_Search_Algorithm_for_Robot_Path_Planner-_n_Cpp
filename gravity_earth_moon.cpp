//Markus Buchholz
//g++ gravity_earth_moon.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>
#include <cmath>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//----------- system dynamic parameters --------------------
const float e = 10;
float G = 6.674 * std::pow(e, -11); // Gravitational constant [m3/(kg s2)]

float M = 5.972 * std::pow(e, 24); // Weight of earth [kg]

float m = 7.348 * std::pow(e, 22); // Weight of moon [kg]

float dt = 10;    // Time step [s]
float v2 = 11200; // Second cosmic velocity

//-----------------------------------------------------------

float radius(float x, float vx, float y, float vy)
{

    return std::pow(std::sqrt(x * x + y * y), 3);
}
//-----------------------------------------------------------
float function1(float x, float vx, float y, float vy)
{

    return vx;
}

//-----------------------------------------------------------
float function2(float x, float vx, float y, float vy)
{

    return -G * M * x / radius(x, vx, y, vy);
}

//-----------------------------------------------------------
float function3(float x, float vx, float y, float vy)
{

    return vy;
}

//-----------------------------------------------------------
float function4(float x, float vx, float y, float vy)
{

    return -G * M * y / radius(x, vx, y, vy);
}

//-----------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> methodRuneKuttaGravity()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;
    std::vector<float> diffEq3;
    std::vector<float> diffEq4;

    std::vector<float> time;

    // init values
    float x1 = 5 * std::pow(e, 6); // x
    float x2 = 0.0;                // v2;                // vx
    float x3 = 0.0;                //  y
    float x4 = 11200;              // v2;                 // vy
    float t = 0.0;                 // init time


    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    diffEq3.push_back(x3);
    diffEq4.push_back(x4);
    time.push_back(t);

    for (int ii = 0; ii < 20000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2, x3, x4);
        float k12 = function2(x1, x2, x3, x4);
        float k13 = function3(x1, x2, x3, x4);
        float k14 = function4(x1, x2, x3, x4);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k23 = function3(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k24 = function4(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k33 = function3(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k34 = function4(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k43 = function3(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k44 = function4(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        x3 = x3 + dt / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);
        x4 = x4 + dt / 6.0 * (k14 + 2 * k24 + 2 * k34 + k44);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        diffEq3.push_back(x3);
        diffEq4.push_back(x4);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, diffEq3, diffEq4, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::tuple<std::vector<float>, std::vector<float>> data1)
{

    std::vector<float> xX1 = std::get<0>(data1);
    std::vector<float> yY1 = std::get<1>(data1);

    plt::plot(xX1, yY1);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::show();
}

//---------------------------------------------------------------

void plot2D2D(std::tuple<std::vector<float>, std::vector<float>> data1, std::tuple<std::vector<float>, std::vector<float>> data2)
{

    int factor = 10000;
    std::vector<float> xX1 = std::get<0>(data1);
    std::vector<float> _yY1 = std::get<1>(data1);

    std::vector<float> yY1;

    for (auto &ii : _yY1)
    {

        yY1.push_back(ii * factor);
    }

    std::vector<float> xX2 = std::get<0>(data2);
    std::vector<float> yY2 = std::get<1>(data2);

    plt::plot(xX1, yY1);
    plt::plot(xX2, yY2);
    plt::xlabel("X[m]");
    plt::ylabel("Y[m]");
    plt::show();
}

//---------------------------------------------------------------
int main()
{
    // x, vx, y, vy, t
    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> dyn = methodRuneKuttaGravity();

    std::tuple<std::vector<float>, std::vector<float>> trajXY = std::make_tuple(std::get<0>(dyn), std::get<2>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> vX = std::make_tuple(std::get<4>(dyn), std::get<1>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> X = std::make_tuple(std::get<4>(dyn), std::get<0>(dyn));
    plot2D(trajXY);
    // plot2D(vX);
    // plot2D(X);
    // plot2D2D(vX, X);
}
