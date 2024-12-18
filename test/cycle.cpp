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
    // auto ghostNeighbour = GhostNeighbour({1, 3}, grid);

    // auto ghostNeighbour2 = GhostNeighbour({12, 14}, grid);
    // ghostNeighbour.insert(ghostNeighbour.end(), ghostNeighbour2.begin(), ghostNeighbour2.end());
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
    MultiIndex index{11,10};
    // index = {5, 13}; //2
    index  = {5, 13};
    CutCellMapping cutCellInfo = cellDivision.getCutCellInfo();
    ParaSet para = cutCellInfo[index];
    double lambda_1 = *para.begin();
    std::cout << "Para 1 : " << lambda_1 << std::endl;
    double lambda_2 = *para.rbegin();
    std::cout << "Para 2 : " << lambda_2 << std::endl;
    ParaInterval lambda = std::make_pair(lambda_1, lambda_2);

    CutCellHandler cutCellHandler;
    cutCellHandler.ClassifyCorners(index, grid, cycle);
    VecList outside_corners = cutCellHandler.getOutsideCorners();
    VecList inside_corners = cutCellHandler.getInsideCorners();
    std::cout << "Inside Corner Size: " << inside_corners.size() << std::endl;
    std::cout << "Outside Corner Size: " << outside_corners.size() << std::endl;

    cutCellHandler.Handler(index, para, grid, cycle);
    double volume = cutCellHandler.getVolume();
    std::cout << "Volume: " << volume << std::endl;
    Vec centroid = cutCellHandler.getCentroid();
    std::cout << "Centroid: " << centroid << std::endl;
    // VecList outside_corners = cutCellHandler.getInsideCorners(index, grid, cycle);
    // std::cout << "Outside Corner Size: " << outside_corners.size() << std::endl;
    // for(const auto& corner : outside_corners)
    // {
    //     std::cout << corner << std::endl;
    // }

    // CurvedTriangle curvedTriangle(cycle, outside_corners[0], lambda);
    // double volume_curved = curvedTriangle.getVolume();
    // std::cout << "Curved Volume: " << volume_curved << std::endl;

    std::cout << "Cell Volume: " << grid.get_cell_volume() << std::endl;
    // std::cout << "Test Handler Volume: " << grid.get_cell_volume() - volume << std::endl;

    // double diff = volume - volume_curved;
    // std::cout << "Diff: " << diff << std::endl;
    // Vec centroid_1 = quad2D_vector(curvedTriangle., lambda_1, lambda_2, 0, 1) / quadresult;
    // std::cout << "Centroid: " << centroid << std::endl;
    // double distance = norm(centroid - cycle.getCenter()) - cycle.getRadius();
    // std::cout << "Distance: " << distance << std::endl;
    plt::axis("equal");
    plotVec(centroid);
    plotGrid(grid);  
    plotCycle(cycle);
    plotCells(cellDivision, grid);
    plt::save("cycle.png");

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