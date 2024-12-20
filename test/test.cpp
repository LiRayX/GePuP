#include "../src/Grid.h"
#include "../src/numlib.h"
#include "../src/Vec.h"
#include "../src/Function.h"
#include "../src/CyclicCurve/CyclicCurve.h"
#include "../src/CyclicCurve/CellDivision.h"
#include "../src/CyclicCurve/CurvedTriangle.h"
#include "../src/CyclicCurve/CutCellHandler.h"
#include "../src/Data.h"

#include <iostream>
class testfun: public ScalarFunction
{
public:
    double operator()(Vec _v) const override
    {
        return 1.0;
    }
};



int main()
{
    Vec low{0, 0};
    Vec high{1, 1};
    int n_seg = 4;
    double h = 1.0 / n_seg;
    Grid grid(low, high, h);
    MultiIndexSet outer_cells = getOuterCells(grid);
    for (const auto &cell : outer_cells)
    {
        std::cout << cell << std::endl;
    }
    std::cout << outer_cells.size() << std::endl;
    std::cout << "----------------" << std::endl;
    MultiIndexList outer_cells_list = grid.getOuterCells();
    for (const auto &cell : outer_cells_list)
    {
        std::cout << cell << std::endl;
    }
    std::cout << outer_cells_list.size() << std::endl;
    // VectorData data(grid);
    // data(6,14) = Vec{0.4375,0.9375};
    // Vec vertex = data(6,14);
    // std::cout << data << std::endl;

    // std::cout << data(6,14) << std::endl;


    //圆
    // Vec center{0.52, 0.51+4*h};
    // double radius = 0.2;
    // CyclicCurve cycle(center, radius);
    // cycle.setIntersections(grid);

    // CellDivision cellDivision;
    // cellDivision.LocateAllCells(grid, cycle);
    // MultiIndex index{6,14};


    // CutCellMapping cutCellInfo = cellDivision.getCutCellInfo();
    // ParaSet para = cutCellInfo[index];
    // Vec vertex_1 = {0.4375,0.9375};
    // Vec In_1 = {0.427838,0.9375};
    // // Vec vertex = {0.4375,0.875};
    // Triangle triangle_1(vertex_1, vertex, In_1);
  
}