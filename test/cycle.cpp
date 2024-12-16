#include "../src/CyclicCurve/CurvedTriangle.h"
#include "../src/CyclicCurve//CellDivision.h"
#include "../src/CyclicCurve/CyclicCurve.h"
#include "../src/CyclicCurve/CutCellHandler.h"

#include "../src/Plot.h"
#include "../src/numlib.h"
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
    auto ghostNeighbour = GhostNeighbour({1, 3}, grid);

    auto ghostNeighbour2 = GhostNeighbour({12, 14}, grid);
    ghostNeighbour.insert(ghostNeighbour.end(), ghostNeighbour2.begin(), ghostNeighbour2.end());
    // for(const auto& index : ghostNeighbour)
    // {
    //     std::cout << index[0] << " " << index[1] << std::endl;
    // }

    // for(const auto& index : ghostNeighbour2)
    // {
    //     std::cout << index[0] << " " << index[1] << std::endl;
    // }
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
    MultiIndex index{6,14};
    CutCellMapping cutCellInfo = cellDivision.getCutCellInfo();
    ParaSet para = cutCellInfo[index];

    CutCellHandler cutCellHandler;
    std::cout << "ParaSet Size: " << para.size() << std::endl;
    double lambda_1 = *para.begin();
    double lambda_2 = *para.rbegin();

    // VecList outside_corners = cutCellHandler.getInsideCorners(index, grid, cycle);
    // std::cout << "Outside Corner Size: " << outside_corners.size() << std::endl;
    // for(const auto& corner : outside_corners)
    // {
    //     std::cout << corner << std::endl;
    // }

    // CurvedTriangle curvedTriangle(cycle, outside_corners[0]);
    // double quadresult = quad2D_scalar(curvedTriangle.Jacobian(), lambda_1, lambda_2, 0, 1); 
    // std::cout << "Quad Result: " << quadresult << std::endl;
    // Vec centroid = quad2D_vector(curvedTriangle.Centroid(), lambda_1, lambda_2, 0, 1) / quadresult;
    // std::cout << "Centroid: " << centroid << std::endl;
    // double distance = norm(centroid - cycle.getCenter()) - cycle.getRadius();
    // std::cout << "Distance: " << distance << std::endl;

    MultiIndex index2{11,10};
    ParaSet para2 = cutCellInfo[index2];
    double lambda_3 = *para2.begin();
    double lambda_4 = *para2.rbegin();
    double lambda[2] = {lambda_3, lambda_4};
    VecList inside_corners2 = cutCellHandler.getInsideCorners(index2, grid, cycle);
    std::cout << "Inside Corner Size: " << inside_corners2.size() << std::endl;
    CurvedTriangle curvedTriangle2(cycle, inside_corners2[0], lambda);
    // std::cout << "Centroid2: " << centroid2 << std::endl;
    // double distance2 = norm(centroid2 - cycle.getCenter()) - cycle.getRadius();
    // std::cout << "Distance2: " << distance2 << std::endl;
    double volume = grid.get_cell_volume() - curvedTriangle2.getVolume();
    std::cout << "Volume: " << volume << std::endl;
    Vec centroid3 = (grid.center(index2)*grid.get_cell_volume() - curvedTriangle2.getCentroid()*curvedTriangle2.getVolume()) / volume;
    std::cout << "Centroid3: " << centroid3 << std::endl;
    double distance3 = norm(centroid3 - cycle.getCenter()) - cycle.getRadius();
    std::cout << "Distance3: " << distance3 << std::endl;

    VecList testEmpty;
    std::cout << "Empty size: " << testEmpty.size() << std::endl;

    // plotVec(centroid);
    plotVec(centroid3);
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