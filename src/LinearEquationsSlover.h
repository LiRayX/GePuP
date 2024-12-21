#pragma once

#include "Vec.h"
#include "Grid.h"
#include "MultiIndexSet.h"
#include "Data.h"



ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, int max_iter);
ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, double tol);


ScalarData Interpolation(const ScalarData &data, const Grid &current_grid, const Grid &refine_grid);
ScalarData Restriction(const ScalarData &data, const Grid &current_grid, const Grid &coarse_grid);


ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, int max_iter)
{
    Eigen::VectorXd D = A.diagonal();
    Eigen::VectorXd D_inv = D.cwiseInverse();
    ScalarData x = x0;
    for(int k=0; k<max_iter; k++)
    {
        Eigen::VectorXd residual = (b - Multiply(A,x)).get_data();
        x.get_data() += weight * D_inv.cwiseProduct(residual);
    }
    return x;
}

ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, double tol)
{
    Eigen::VectorXd D = A.diagonal();
    Eigen::VectorXd D_inv = D.cwiseInverse();
    ScalarData x = x0;
    Eigen::VectorXd residual = (b - Multiply(A,x)).get_data();
    while(residual.norm() > tol)
    {
        x.get_data() += weight * D_inv.cwiseProduct(residual);
        residual = (b - Multiply(A,x)).get_data();
    }
    return x;
}


ScalarData Interpolation(const ScalarData &data, const Grid &current_grid, const Grid &refine_grid)
{
    ScalarData result = ScalarData(refine_grid);
    loop_cell_2(current_grid, i, j)
    {
        MultiIndex current_cell{i, j};
        MultiIndex index{2*i, 2*j};
        MultiIndexList refine_index = {index, index + MultiIndex{1, 0}, index + MultiIndex{0, 1}, index + MultiIndex{1, 1}};
        for(const auto &refine_cell : refine_index)
        {
            result(refine_cell) = data(current_cell);
        }
    }
    return result;
}

ScalarData Restriction(const ScalarData &data, const Grid &current_grid, const Grid &coarse_grid)
{
    ScalarData result = ScalarData(coarse_grid);
    loop_cell_2(coarse_grid, i, j)
    {
        MultiIndex coarse_cell{i, j};
        MultiIndex index{2*i, 2*j};
        MultiIndexList current_index = {index, index + MultiIndex{1, 0}, index + MultiIndex{0, 1}, index + MultiIndex{1, 1}};
        double sum = 0;
        for(const auto &current_cell : current_index)
        {
            sum += data(current_cell);
        }
        result(coarse_cell) = sum / 4;
    }
    return result;
}


ScalarBoundaryFaceAvr FaceAvrRestriction(const ScalarBoundaryFaceAvr &data, const Grid &current_grid, const Grid &coarse_grid)
{
    ScalarBoundaryFaceAvr result = ScalarBoundaryFaceAvr(coarse_grid);
    MultiIndexSet boundaryCells = getBoundaryCells(coarse_grid);
    for(int k=0; k<coarse_grid.get_size()[1];k++)
    {
        result.get_face_avr_l()(k) = data.get_face_avr_l()(2*k) + data.get_face_avr_l()(2*k+1);
        result.get_face_avr_r()(k) = data.get_face_avr_r()(2*k) + data.get_face_avr_r()(2*k+1);
    }
    for(int k=0; k<coarse_grid.get_size()[0];k++)
    {
        result.get_face_avr_d()(k) = data.get_face_avr_d()(2*k) + data.get_face_avr_d()(2*k+1);
        result.get_face_avr_u()(k) = data.get_face_avr_u()(2*k) + data.get_face_avr_u()(2*k+1);
    }
    return result;
}



// ScalarData Vcycle(const Grid &grid, const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, int level, int weight, int max_iter)
// {
//     Grid current_grid = grid;
//     if(level == 0)
//     {
//         return WeightedJacobi(A, b, x0, weight, max_iter);
//     }
//     else
//     {
//         ScalarData x = WeightedJacobi(A, b, x0, weight, max_iter);
//         ScalarData residual = b - Multiply(A, x);
//         Grid grid = current_grid.coarsen();
//         ScalarData error = Vcycle(grid, A, residual, ScalarData(grid), level-1, weight, max_iter);
//         x += error;
//         return x;
//     }
// }
