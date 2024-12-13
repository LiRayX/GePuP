// Draw a simple case using matplotlib-cpp
#include "BoundaryCurve.h"
#include "Grid.h"
#include "Vec.h"
#include "../lib/matplotlib-cpp/matplotlibcpp.h"
#include "CurveBelonging.h"
#include "CyclicCurve.h"
#include "CellDivision.h"
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
void plotPieceWiseSpline(const Spline2d &spline, const CurveBelonging &curveBelonging, double step);
void plotWholeBoundaryCurve(const BoundaryCurve &curve);

void plotCycle(const CyclicCurve &cycle);
void plotCells(const CellDivision &cellDivison, const Grid &grid);
void plotCellsWithLegend(const CellDivision &cellDivison, const Grid &grid);



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
        // plt::plot(x_grid, y_grid);
        plt::plot(x_grid, y_grid, "k-");
}

void plotControlPoints(const VecList &controlPoints)
{
    std::vector<double> x, y;
    for (const auto &point : controlPoints)
    {
        x.push_back(point[0]);
        y.push_back(point[1]);
    }
    plt::scatter(x, y, 5.0, {{"color", "red"}});
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

void plotPieceWiseSpline(const Spline2d &spline, const CurveBelonging &curveBelonging, double step)
{
    const ParaIntervalList &paraIntervals = curveBelonging.getParaIntervals();
    for (size_t i = 0; i < paraIntervals.size(); ++i)
    {
        double start = paraIntervals[i].first;
        double end = paraIntervals[i].second;
        plotSpline(spline, start, end, step);
    }
}

void plotIntersectionPoints(const Spline2d &spline, const CurveBelonging &curveBelonging)
{
    std::vector<double> x, y;
    VecList intersection_list = curveBelonging.getIntersectionPoints(spline);
    for (const auto &point : intersection_list)
    {
        x.push_back(point[0]);
        y.push_back(point[1]);
    }
    plt::scatter(x, y, 10.0, {{"color", "purple"}, {"marker", "x"}});
}

void plotWholeBoundaryCurvePieceWisely(const BoundaryCurve &curve)
{   //plot control points
    const VecList &controlPoints = curve.getControlPoints();
    plotControlPoints(controlPoints);
    //plot splines and intersection points
    const SplineList &splines = curve.getSplines();
    const PieceWiseBelongingList &PieceWiseBelongingIfo = curve.getPieceWiseBelongingIfo();
    for(size_t i = 0; i < splines.size(); i++)
    {
        plotPieceWiseSpline(splines[i], PieceWiseBelongingIfo[i], 0.01);
        plotIntersectionPoints(splines[i], PieceWiseBelongingIfo[i]);
    }
}

void plotWholeBoundaryCurve(const BoundaryCurve &curve)
{   //plot control points
    const VecList &controlPoints = curve.getControlPoints();
    plotControlPoints(controlPoints);
    //plot splines
    const SplineList &splines = curve.getSplines();
    plotSplineList(splines, 0, 1, 0.01);
    //plot intersection points
    const PieceWiseBelongingList &PieceWiseBelongingIfo = curve.getPieceWiseBelongingIfo();
    for(size_t i = 0; i < splines.size(); i++)
    {
        plotIntersectionPoints(splines[i], PieceWiseBelongingIfo[i]);
    }
}

void plotCycle(const CyclicCurve &cycle)
{
    Vec center = cycle.getCenter();
    double radius = cycle.getRadius();
    double theta = 0;
    std::vector<double> x, y;
    for (theta = 0; theta <= 2 * Pi; theta += 0.01)
    {
        x.push_back(center[0] + radius * std::cos(theta));
        y.push_back(center[1] + radius * std::sin(theta));
    }
    x.push_back(center[0] + radius);
    y.push_back(center[1]);
    plt::plot(x, y, "b-");
}

void plotCellsWithLegend(const CellDivision &cellDivison, const Grid &grid)
{
    std::vector<double> centerX, centerY;

    // 绘制 DeadCells
    const MultiIndexSet &DeadCells = cellDivison.getDeadCells();
    for (const auto &index : DeadCells)
    {
        Vec center = grid.center(index[0], index[1]);
        centerX.push_back(center[0]);
        centerY.push_back(center[1]);
    }
    plt::named_plot("DeadCells", centerX, centerY, "ko");
    centerX.clear();
    centerY.clear();

    // 绘制 CutCells
    const MultiIndexSet &CutCells = cellDivison.getCutCells();
    for (const auto &index : CutCells)
    {
        Vec center = grid.center(index[0], index[1]);
        centerX.push_back(center[0]);
        centerY.push_back(center[1]);
    }
    plt::named_plot("CutCells", centerX, centerY, "bo");
    centerX.clear();
    centerY.clear();

    // 绘制 SideCells
    const MultiIndexSet &SideCells = cellDivison.getSideCells();
    for (const auto &index : SideCells)
    {
        Vec center = grid.center(index[0], index[1]);
        centerX.push_back(center[0]);
        centerY.push_back(center[1]);
    }
    plt::named_plot("SideCells", centerX, centerY, "go");
    centerX.clear();
    centerY.clear();

    // 绘制 EdgeCells
    const MultiIndexSet &EdgeCells = cellDivison.getEdgeCells();
    for (const auto &index : EdgeCells)
    {
        Vec center = grid.center(index[0], index[1]);
        centerX.push_back(center[0]);
        centerY.push_back(center[1]);
    }
    plt::named_plot("EdgeCells", centerX, centerY, "yo");
    centerX.clear();
    centerY.clear();

    // 绘制 CoreCells
    const MultiIndexSet &CoreCells = cellDivison.getCoreCells();
    for (const auto &index : CoreCells)
    {
        Vec center = grid.center(index[0], index[1]);
        centerX.push_back(center[0]);
        centerY.push_back(center[1]);
    }
    plt::named_plot("CoreCells", centerX, centerY, "mo");

    // 显示图例
    plt::legend();
}


void plotCells(const CellDivision &cellDivison, const Grid &grid)
{
    //plot DeadCells
    const MultiIndexSet &DeadCells = cellDivison.getDeadCells();
    for (const auto &index : DeadCells)
    {
        Vec center = grid.center(index[0], index[1]);
        plt::plot({center[0]}, {center[1]}, "ko");
    }
    //plot AliveCells
    //plot CutCells
    const MultiIndexSet &CutCells = cellDivison.getCutCells();
    for (const auto &index : CutCells)
    {
        Vec center = grid.center(index[0], index[1]);
        plt::plot({center[0]}, {center[1]}, "bo");
    }
    //plot SideCells
    const MultiIndexSet &SideCells = cellDivison.getSideCells();
    for (const auto &index : SideCells)
    {
        Vec center = grid.center(index[0], index[1]);
        plt::plot({center[0]}, {center[1]}, "go");
    }
    //plot EdgeCells
    const MultiIndexSet &EdgeCells = cellDivison.getEdgeCells();
    for (const auto &index : EdgeCells)
    {
        Vec center = grid.center(index[0], index[1]);
        plt::plot({center[0]}, {center[1]}, "yo");
    }
    // plot CoreCells
    const MultiIndexSet &CoreCells = cellDivison.getCoreCells();
    for (const auto &index : CoreCells)
    {
        Vec center = grid.center(index[0], index[1]);
        plt::plot({center[0]}, {center[1]}, "mo");
    }
}
