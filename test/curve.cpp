// Draw a simple case using matplotlib-cpp
#include "../src/BoundaryCurve.h"
#include "../src/Grid.h"
#include "../src/Vec.h"
#include "../src/CurveBelonging.h"
#include "../src/Plot.h"
#include "../src/CellClassifier.h"
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
    VecList endPointsVec;
    for(const auto& index : endPoints)
    {
        endPointsVec.push_back(controlPoints[index]);
    }
    BoundaryCurve curve(controlPoints, endPoints);
    SplineList splines = curve.getSplines();
    // Spline2d spline  = splines[0];
    // const Spline2d &spline = curve.getSplines()[0];

    double tol = 1e-6;
    double para_tol = 1e-12;
    int max_iter = 100;
    curve.setPieceWiseBelongingIfo(g, tol, para_tol, max_iter);
    // const PieceWiseBelongingList &PieceWiseBelongingIfo = curve.getPieceWiseBelongingIfo();
    // const CurveBelonging &curveBelonging = PieceWiseBelongingIfo[0];
    // int n_cutcell = curveBelonging.getMultiIndices().size();
    // std::cout << "MultiIndices Size: " << n_cutcell << std::endl;
    CellClassifier cellClassifier;
    cellClassifier.LocateCutCells(g, curve);
    const MultiIndexSet& cutCells = cellClassifier.getCutCells();
    const CutCellMap& cutCellInfo = cellClassifier.getCutCellInfo();
    std::cout << "CutCells: " << std::endl;
    for(const auto& index : cutCells)
    {
        std::cout << index[0] << " " << index[1] << std::endl;
    }
    std::cout << "CutCellInfo: " << std::endl;
    for(const auto& info : cutCellInfo)
    {
        std::cout << "Index: " << info.first[0] << " " << info.first[1] << std::endl;
        std::cout << "IndexSpline: " << std::endl;
        for(const auto& index : info.second.index_spline)
        {
            std::cout << index << std::endl;
        }
        std::cout << "ParaInterval: " << std::endl;
        for(const auto& interval : info.second.parainterval)
        {
            std::cout << "[" << interval.first << " ," << interval.second << "]" << std::endl;
        }
    }




    /****************************************************************************/



    // 绘制网格线
    plotGrid(g);
    // 绘制边界曲线
    plotWholeBoundaryCurve(curve);
    // // 绘制控制点
    // plotControlPoints(controlPoints);
    // //绘制曲线
    // plotPieceWiseSpline(spline, curveBelonging, 0.01);
    // //绘制交点
    // plotIntersectionPoints(spline, curveBelonging);
    //显示图像
    plt::save("whole_curve.svg");







    return 0;


    
}