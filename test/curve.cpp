// Draw a simple case using matplotlib-cpp
#include "../src/BoundaryCurve.h"
#include "../src/Grid.h"
#include "../src/Vec.h"
#include "../src/CurveBelonging.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>


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
    Spline2d spline  = splines[0];

    CurveBelonging curveBelonging;
    curveBelonging.Bisection(g, spline, 0.05, 1e-6);
    const MultiIndexList& multiIndices = curveBelonging.getMultiIndices();
    const ParaIntervalList& paraIntervals = curveBelonging.getParaIntervals();
    const MultiIndexSet& localCutCells = curveBelonging.getLocalCutCells();
    std::cout << "MultiIndices: " << std::endl;
    for(const auto& index : multiIndices)
    {
        std::cout << index[0] << " " << index[1] << std::endl;
    }
    std::cout << "ParaIntervals: " << std::endl;
    for(const auto& interval : paraIntervals)
    {
        std::cout << interval.first << " " << interval.second << std::endl;
    }
    std::cout << "LocalCutCells: " << std::endl;
    for(const auto& index : localCutCells)
    {
        std::cout << index[0] << " " << index[1] << std::endl;
    }
    return 0;
}