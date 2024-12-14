#include "../src/CyclicCurve.h"
#include "../src/CellDivision.h"
#include "../src/Plot.h"

#include <iostream>



int main()
{
     //网格信息
    Vec low{0, 0};
    Vec high{1, 1};
    int n_seg = 16;
    double h = 1.0 / n_seg;
    Grid grid(low, high, h);
    //圆
    Vec center{0.52, 0.51+4*h};
    double radius = 0.2;
    CyclicCurve cycle(center, radius);
    cycle.setIntersections(grid);

    CellDivision cellDivision;
    cellDivision.LocateAllCells(grid, cycle);

    const MultiIndexSet& cutCells = cellDivision.getCutCells();
    const MultiIndexSet& sideCells = cellDivision.getSideCells();
    // std::cout << "sideCells Size: " << sideCells.size() << std::endl;
    // for (const auto& index : sideCells)
    // {
    //     std::cout << index[0] << " " << index[1] << std::endl;
    // }
    // std::cout << "CutCells Size: " << cutCells.size() << std::endl;


    // for(const auto& index : cutCells)
    // {
    //     std::cout << index[0] << " " << index[1] << std::endl;
    // }
    // auto cutCellInfo = cellDivision.getCutCellInfo();
    // for(const auto& info : cutCellInfo)
    // {
    //     std::cout << "Index: " << info.first[0] << " " << info.first[1] << std::endl;
    //     for(const auto& para : info.second)
    //     {
    //         std::cout << "Para: " << para << std::endl;
    //     }
    // }

    plotGrid(grid);  
    plotCycle(cycle);
    plotCells(cellDivision, grid);
    plt::save("cycle.svg");

    // std::cout << "Intersections Size: " << intersections.size() << std::endl;
    // for(const auto& para : intersections)
    // {
    //     std::cout << "Para: " << para << std::endl;
    // }

    // std::cout << "Intersections Size: " << intersections2.size() << std::endl;
    // for(const auto& para : intersections2)
    // {
    //     std::cout << "Para: " << para << std::endl;
    // }










    // Vec point{1.0, 0.0};
    // Vec normal = curve.getNormal(point);
    // Vec tangent = curve.getTangent(point);
    // std::cout << "Normal: " << normal[0] << " " << normal[1] << std::endl;
    // std::cout << "Tangent: " << tangent[0] << " " << tangent[1] << std::endl;

    // double angle = Pi/4;
    // Vec point2 = curve.getPoint(angle);
    // Vec normal2 = curve.getNormal(angle);
    // Vec tangent2 = curve.getTangent(angle);
    // std::cout << "Point: " << point2[0] << " " << point2[1] << std::endl;
    // std::cout << "Normal: " << normal2[0] << " " << normal2[1] << std::endl;
    // std::cout << "Tangent: " << tangent2[0] << " " << tangent2[1] << std::endl;

    // std::cout<<"Test Arccos"<< std::acos(0.707107) << std::endl;
    // std::cout<<"Test Arccos"<< std::acos(-0.707107) << std::endl;
    // std::cout<<"Test Arcsin"<< std::asin(0.707107) << std::endl;
    // std::cout<<"Test Arcsin"<< std::asin(-0.707107) << std::endl;
}