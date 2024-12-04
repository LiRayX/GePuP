// Draw a simple case using matplotlib-cpp
#include "../src/BoundaryCurve.h"
#include "../src/Grid.h"
#include "../src/Vec.h"
#include "../src/CurveBelonging.h"
#include "../src/Plot.h"
#include "../lib/matplotlib-cpp/matplotlibcpp.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>


int main()
{
     //网格信息
    Vec low{0, 0};
    Vec high{1, 1};
    int n_seg = 8;
    double h = 1.0 / n_seg;
    Grid g(low, high, h);
    //控制点List
    VecList controlPoints;
    //控制点按照逆时针方向
    int num_points = 30;
    double radius = 0.26;
    for(int i = 0; i <= num_points; i++)
    {
        double x = 0.5 + radius * cos(2 * M_PI * i / num_points);
        double y = 0.5 + radius * sin(2 * M_PI * i / num_points);
        controlPoints.push_back(Vec{x, y});
    }
    //边界曲线
    std::vector<int> endPoints;
    for(int i = 0; i <= num_points; i+=3)
    {
        endPoints.push_back(i);
    }
    BoundaryCurve curve(controlPoints, endPoints);
    SplineList splines = curve.GetSplines();
    Spline2d spline  = splines[1];
    CurveBelonging curveBelonging;
    curveBelonging.AdaptiveCheck(g, spline);
    
    const MultiIndexList& multiIndices = curveBelonging.getMultiIndices();
    const ParaIntervalList& paraIntervals = curveBelonging.getParaIntervals();
    const MultiIndexSet& localCutCells = curveBelonging.getLocalCutCells();
    // Normal normal = minus(multiIndices[1] , multiIndices[0]);
    double tol = 1e-6;
    double para_tol = 1e-20;
    int max_iter = 100;
    // double bisection_result = IntervalBisection(g, spline, paraIntervals[0], multiIndices[0], Normal{1, 0}, tol, para_tol, 100);
    // std::cout << "First Bisection Result: " << bisection_result << std::endl;
    // Vec intersection{spline(bisection_result)};
    // std::cout << "Current Cell: " << std::endl;
    // std::cout << multiIndices[0][0] << " " << multiIndices[0][1] << std::endl;
    // std::cout << "Next Cell: " << std::endl;
    // std::cout << multiIndices[1][0] << " " << multiIndices[1][1] << std::endl;
    // std::cout << "Normal: " << std::endl;
    // std::cout << normal[0] << " " << normal[1] << std::endl;
    // std::cout << "test of Interval Bisection " << std::endl;
    // std::cout << "Bisection Result: " << bisection_result << std::endl;
    // std::cout << "Intersection: " << intersection << std::endl;
    VecList intersection_list = curveBelonging.getIntersectionPoints(spline);
    std::cout << "Before Bisection " << std::endl;
    std::cout << "MultiIndices: " << std::endl;
    for(const auto& index : multiIndices)
    {
        std::cout << index[0] << " " << index[1] << std::endl;
    }
    std::cout << "ParaIntervals: " << std::endl;
    for(const auto& interval : paraIntervals)
    {
        std::cout <<"[" << interval.first << " ," << interval.second << "]" << std::endl;
    } 
    std::cout << "Intersection Points: " << std::endl;
    for(const auto& intersection : intersection_list)
    {
        std::cout << intersection << std::endl;
    }
    std::cout << "LocalCutCells: " << std::endl;
    for(const auto& index : localCutCells)
    {
        std::cout << index[0] << " " << index[1] << std::endl;
    }

    curveBelonging.PieceWiseBelonging(g, spline, tol, para_tol, max_iter);
    intersection_list.clear();
    intersection_list = curveBelonging.getIntersectionPoints(spline);
    std::cout << "After Bisection " << std::endl;
    std::cout << "ParaIntervals: " << std::endl;
    for(const auto& interval : paraIntervals)
    {
        std::cout <<"[" << interval.first << " ," << interval.second << "]" << std::endl;
    } 
    std::cout << "Intersection Points: " << std::endl;
    for(const auto& intersection : intersection_list)
    {
        std::cout << intersection << std::endl;
    }
    // 绘制网格线
    plotGrid(g);
    // 绘制控制点
    plotControlPoints(controlPoints);
    //绘制曲线
    plotPieceWiseSpline(spline, curveBelonging, 0.01);
    //绘制交点
    plotIntersectionPoints(spline, curveBelonging);
    //显示图像
    plt::save("piecewise_spline.svg");







    return 0;


    
}