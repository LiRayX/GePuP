// Draw a simple case using matplotlib-cpp
#include "../src/BoundaryCurve.h"
#include "../src/Grid.h"
#include "../src/Vec.h"
// #include "../lib/Eigen/Core"
// #include "../lib/unsupported/Eigen/Splines"
// #include "../lib/Eigen/Dense"
#include "../lib/matplotlib-cpp/matplotlibcpp.h"
#include <iostream>
#include <vector>
#include <cmath>

namespace plt = matplotlibcpp;

int main()
{
    //网格信息
    Vec low{0, 0};
    Vec high{1, 1};
    int n_seg = 4;
    double h = 1.0 / n_seg;
    Grid g(low, high, h);
    //网格线
    std::vector<double> x_grid, y_grid;
    for (int i = 0; i <= n_seg; ++i)
    {
        double coord = i * h;
        // 垂直线
        x_grid.push_back(coord);
        y_grid.push_back(0);
        x_grid.push_back(coord);
        y_grid.push_back(1);
        x_grid.push_back(std::numeric_limits<double>::quiet_NaN()); // 分隔线段
        y_grid.push_back(std::numeric_limits<double>::quiet_NaN()); // 分隔线段

        // 水平线
        x_grid.push_back(0);
        y_grid.push_back(coord);
        x_grid.push_back(1);
        y_grid.push_back(coord);
        x_grid.push_back(std::numeric_limits<double>::quiet_NaN()); // 分隔线段
        y_grid.push_back(std::numeric_limits<double>::quiet_NaN()); // 分隔线段
    }
    //控制点List
    VecList controlPoints;
    //控制点按照逆时针方向
    int num_points = 10;
    double radius = 0.25;
    for(int i = 0; i <= num_points; i++)
    {
        double x = 0.5 + radius * cos(2 * M_PI * i / num_points);
        double y = 0.5 + radius * sin(2 * M_PI * i / num_points);
        controlPoints.push_back(Vec{x, y});
    }


    // 绘制网格线
    plt::plot(x_grid, y_grid);

    // 绘制控制点
    std::vector<double> control_x, control_y;
    for (const auto &point : controlPoints)
    {
        control_x.push_back(point[0]);
        control_y.push_back(point[1]);
    }
    plt::scatter(control_x, control_y, 10.0, {{"color", "red"}});



    plt::save("output.png");
    return 0;
}