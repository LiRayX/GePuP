// Draw a simple case using matplotlib-cpp
#include "../src/BoundaryCurve.h"
#include "../src/Grid.h"
#include "../src/Vec.h"
#include "../src/Plot.h"
// #include "../lib/Eigen/Core"
// #include "../lib/unsupported/Eigen/Splines"
#include "../lib/matplotlib-cpp/matplotlibcpp.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

namespace plt = matplotlibcpp;

// /// @brief plot grid, control points and splines.
// void plotGrid(const Grid &g);
// void plotControlPoints(const VecList &controlPoints);
// void plotSpline(const Spline2d &spline, double start, double end, double step);
// void plotSplineList(const SplineList &splines, double start, double end, double step);

int main()
{
    //网格信息
    Vec low{0, 0};
    Vec high{1, 1};
    int n_seg = 4;
    double h = 1.0 / n_seg;
    Grid g(low, high, h);
    //控制点List
    VecList controlPoints;
    //控制点按照逆时针方向
    int num_points = 9;
    double radius = 0.25;
    for(int i = 0; i <= num_points; i++)
    {
        double x = 0.5 + radius * cos(2 * M_PI * i / num_points);
        double y = 0.5 + radius * sin(2 * M_PI * i / num_points);
        controlPoints.push_back(Vec{x, y});
    }
    //边界曲线
    std::vector<int> endPoints{0, 3, 6, 9};
    BoundaryCurve curve(controlPoints, endPoints);
    SplineList splines = curve.GetSplines();

    // 绘制网格线
    plotGrid(g);
    // 绘制控制点
    plotControlPoints(controlPoints);
    //绘制曲线
    plotSplineList(splines, 0, 1, 0.01);
    //显示图像
    plt::save("output.png");
    return 0;
}

