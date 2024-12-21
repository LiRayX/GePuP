#include "../src/Grid.h"
#include "../src/numlib.h"
#include "../src/Vec.h"
#include "../src/Function.h"
#include "../src/CyclicCurve/CyclicCurve.h"
#include "../src/CyclicCurve/CellDivision.h"
#include "../src/CyclicCurve/CurvedTriangle.h"
#include "../src/CyclicCurve/CutCellHandler.h"
#include "../src/Data.h"
#include "../src/DiffMatrix.h"
#include "../src/LinearEquationsSlover.h"

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
    int n_seg = 2;
    double h = 1.0 / n_seg;
    Grid grid(low, high, h);
    // ScalarBoundaryFaceAvr scalarBoundaryFaceAvr(grid);
    // VectorBoundaryFaceAvr vectorBoundaryFaceAvr(grid);
    // type boundary_type = type::Dirichlet;
    // // LaplacianOperator laplacian(grid, scalarBoundaryFaceAvr, boundary_type);
    // GradientOperator gradient(grid, scalarBoundaryFaceAvr, boundary_type);
    // MatrixXd matrix_x = gradient.get_matrix_x();
    // MatrixXd matrix_y = gradient.get_matrix_y();
    MatrixXd A(4, 4);
    A << 4, -1, 0, 0,
         -1, 4, -1, 0,
         0, -1, 4, -1,
         0, 0, -1, 3;
    VectorXd b(4);
    b << 15, 10, 10, 10;
    ScalarData rhs = ScalarData(grid);
    rhs.set_data(b);
    ScalarData x0 = ScalarData(grid);
    Eigen::SparseMatrix<double> A_S = A.sparseView();
    // ScalarData x = Multiply(A_S, rhs);
    ScalarData x = WeightedJacobi(A_S, rhs, x0, 0.5, 100);
    std::cout << x.get_data() << std::endl;
    // std::cout << matrix_x << std::endl;
    // std::cout<<"-------------------"<<std::endl;
    // std::cout << matrix_y << std::endl;
    // std::cout<<"-------------------"<<std::endl;
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