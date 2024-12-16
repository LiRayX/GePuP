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
    
    CurvedTriangle curvedTriangle(cycle, center, lambda_1, lambda_2);
 /*************************************************************************** */
  //  test for the Jacobian of the CurvedTriangle
  //  Diff: 8.67362e-19
    double area = 0.5* (lambda_2 - lambda_1) * radius * radius;
    std::cout << "Exat Area: " << area << std::endl;
    double quadresult = curvedTriangle.getVolume();
    std::cout << "Quad Result: " << quadresult << std::endl;

    double diff = std::fabs(quadresult - area);
    std::cout << "Diff: " << diff << std::endl;
  
}