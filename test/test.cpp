#include "../src/Grid.h"
#include "../src/numlib.h"
#include "../src/Vec.h"
#include "../src/Function.h"
#include "../src/CyclicCurve/CyclicCurve.h"
#include "../src/CyclicCurve/CellDivision.h"
#include "../src/CyclicCurve/CurvedTriangle.h"
#include "../src/CyclicCurve/CutCellHandler.h"

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
    MultiIndex index{6,14};
    CutCellMapping cutCellInfo = cellDivision.getCutCellInfo();
    ParaSet para = cutCellInfo[index];

    CutCellHandler cutCellHandler;
    std::cout << "ParaSet Size: " << para.size() << std::endl;
    double lambda_1 = *para.begin();
    double lambda_2 = *para.rbegin();

    VecList outside_corners = cutCellHandler.getOutsideCorners(index, grid, cycle);
    std::cout << "Outside Corner Size: " << outside_corners.size() << std::endl;
    for(const auto& corner : outside_corners)
    {
        std::cout << corner << std::endl;
    }

    CurvedTriangle curvedTriangle(cycle, outside_corners[0]);
    double quadresult = quad2D(curvedTriangle.Jacobian(), lambda_1, lambda_2, 0, 1); 
    std::cout << "Quad Result: " << quadresult << std::endl;



  
}