// Draw a simple case using matplotlib-cpp
#include "../src/BoundaryCurve.h"
#include "../src/Grid.h"
#include "../src/Vec.h"
#include "../lib/matplotlib-cpp/matplotlibcpp.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

namespace plt = matplotlibcpp;

/// @brief plot grid, control points and splines.
void plotGrid(const Grid &g);
void plotControlPoints(const VecList &controlPoints);
void plotSpline(const Spline2d &spline, double start, double end, double step);
void plotSplineList(const SplineList &splines, double start, double end, double step);


void plotGrid(const Grid &g)
{
    std::vector<double> x_grid, y_grid;
    assert(g.get_size()[0] == g.get_size()[1]);
    int n_seg = g.get_size()[0];
    double h = g.get_h();
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
    plt::plot(x_grid, y_grid);
}

void plotControlPoints(const VecList &controlPoints)
{
    std::vector<double> x, y;
    for (const auto &point : controlPoints)
    {
        x.push_back(point[0]);
        y.push_back(point[1]);
    }
    plt::scatter(x, y, 10.0, {{"color", "red"}});
}

void plotSpline(const Spline2d &spline, double start, double end, double step)
{
    assert(start <= end && step > 0 && start >= 0 && end <= 1);
    std::vector<double> x, y;
    for (double u = start; u <= end; u += step)
    {
        Eigen::Vector2d point = spline(u);
        x.push_back(point[0]);
        y.push_back(point[1]);
    }
    plt::plot(x, y,"--");
}

void plotSplineList(const SplineList &splines, double start, double end, double step)
{
    for (const auto &spline : splines)
    {
        plotSpline(spline, start, end, step);
    }
}