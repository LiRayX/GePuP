#pragma once

#include "Vec.h"
#include "Grid.h"
#include "MultiIndexSet.h"
#include "Data.h"
#include "Coefficients.h"
#include "../lib/Eigen/Sparse"

using DiffMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
using TripletList = std::vector<Triplet>;
using MultiIndexList = std::vector<MultiIndex>;

class LaplacianOperator
{
public:
    LaplacianOperator(const Grid &grid, const ScalarBoundaryFaceAvr &boudary_value, type boundary_type);
    const DiffMatrix &get_matrix() const { return matrix; }
    const ScalarData &get_rhs() const { return rhs; }
    ~LaplacianOperator() = default;
protected:
    DiffMatrix matrix;
    ScalarData rhs;
};
class GradientOperator
{
public:
    GradientOperator(const Grid &grid, const ScalarBoundaryFaceAvr &boudary_value, type boundary_type);
    const DiffMatrix &get_matrix_x() const { return matrix_x; }
    const DiffMatrix &get_matrix_y() const { return matrix_y; }
    const ScalarData &get_rhs_x() const { return rhs_x; }
    const ScalarData &get_rhs_y() const { return rhs_y; }
    ~GradientOperator() = default;
protected:
    DiffMatrix matrix_x;
    DiffMatrix matrix_y;
    ScalarData rhs_x;
    ScalarData rhs_y;
};
/** 
class DivergenceOperator
{
public:
    DivergenceOperator(const Grid &grid, const VectorBoundaryFaceAvr &boudary_value, type boundary_type);
    const DiffMatrix &get_matrix() const { return matrix; }
    const ScalarData &get_rhs() const { return rhs; }
    ~DivergenceOperator() = default;
protected:
    DiffMatrix matrix;
    ScalarData rhs;
};
*/



GradientOperator::GradientOperator(const Grid &grid, const ScalarBoundaryFaceAvr &boudary_value, type boundary_type): rhs_x(grid), rhs_y(grid) 
{
    int n_cells = grid.get_size()[0] * grid.get_size()[1];
    double h = grid.get_h();
    TripletList tripletList_X;
    TripletList tripletList_Y;
    matrix_x.resize(n_cells, n_cells);
    matrix_y.resize(n_cells, n_cells);

    loop_inner_cell_2(grid,i,j)
    {
        MultiIndex centered_cell{i,j};
        MultiIndexList index_diff_cells_x = FourthOrderDiffCells_X(centered_cell);
        MultiIndexList index_diff_cells_y = FourthOrderDiffCells_Y(centered_cell);
        for (const auto &index : index_diff_cells_x)
        {
            int ManhanDistance = SignManhattanDistance(centered_cell, index);
            double coeff = GradientCoefficients(ManhanDistance,h);
            tripletList_X.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
        }
        for (const auto &index : index_diff_cells_y)
        {
            int ManhanDistance = SignManhattanDistance(centered_cell, index);
            double coeff = GradientCoefficients(ManhanDistance,h);
            tripletList_Y.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
        }
    }
    MultiIndexSet outer_cells = getOuterCells(grid);
    for (const auto &centered_cell : outer_cells)
    {
        MultiIndexList index_diff_cells_x = FourthOrderDiffCells_X(centered_cell);
        MultiIndexList index_diff_cells_y = FourthOrderDiffCells_Y(centered_cell);
        for (const auto &index : index_diff_cells_x)
        {
            if (grid.isIndexValid(index))
            {
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double coeff = GradientCoefficients(ManhanDistance,h);
                tripletList_X.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
            }
            //ghost cell
            else
            {
                //get the coffe of DiffOperator
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double diff_coeff = GradientCoefficients(ManhanDistance,h);
                //get the coeff of GhostCell
                int layer = grid.layer(index);
                Normal inner_normal = centered_cell - index;
                Normal direction = normalize(inner_normal);
                //start from the boundary cell
                MultiIndex boundary_cell = index + direction*layer;
                for(int k=0; k<3; k++)
                {
                    double ghost_coeff = GhostCellCoefficients(layer, k, h, boundary_type);
                    tripletList_X.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(boundary_cell+direction*k), ghost_coeff*diff_coeff));
                }
                rhs_x(centered_cell) += diff_coeff*GhostCellCoefficients(layer, -1, h, boundary_type)*boudary_value(boundary_cell, -1*direction);
            }
        }
        for (const auto &index : index_diff_cells_y)
        {
            if (grid.isIndexValid(index))
            {
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double coeff = GradientCoefficients(ManhanDistance,h);
                tripletList_Y.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
            }
            //ghost cell
            else
            {
                //get the coffe of DiffOperator
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double diff_coeff = GradientCoefficients(ManhanDistance,h);
                //get the coeff of GhostCell
                int layer = grid.layer(index);
                Normal inner_normal = centered_cell - index;
                Normal direction = normalize(inner_normal);
                //start from the boundary cell
                MultiIndex boundary_cell = index + direction*layer;
                for(int k=0; k<3; k++)
                {
                    double ghost_coeff = GhostCellCoefficients(layer, k, h, boundary_type);
                    tripletList_Y.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(boundary_cell+direction*k), ghost_coeff*diff_coeff));
                }
                rhs_y(centered_cell) += diff_coeff*GhostCellCoefficients(layer, -1, h, boundary_type)*boudary_value(boundary_cell, -1*direction);
            }
        }
    }
    matrix_x.setFromTriplets(tripletList_X.begin(), tripletList_X.end());
    matrix_y.setFromTriplets(tripletList_Y.begin(), tripletList_Y.end());
}












LaplacianOperator::LaplacianOperator(const Grid &grid, const ScalarBoundaryFaceAvr &boudary_value, type boundary_type): rhs(grid) 
{
    int n_cells = grid.get_size()[0] * grid.get_size()[1];
    double h = grid.get_h();
    TripletList tripletList;
    matrix.resize(n_cells, n_cells);

    loop_inner_cell_2(grid,i,j)
    {
        MultiIndex centered_cell{i,j};
        MultiIndexList index_diff_cells = FourthOrderDiffCells(centered_cell);
        for (const auto &index : index_diff_cells)
        {
            int ManhanDistance = SignManhattanDistance(centered_cell, index);
            double coeff = LaplacainCoefficients(ManhanDistance,h);
            tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
        }
    }
    MultiIndexSet outer_cells = getOuterCells(grid);
    for (const auto &centered_cell : outer_cells)
    {
        MultiIndexList index_diff_cells = FourthOrderDiffCells(centered_cell);
        for (const auto &index : index_diff_cells)
        {
            if (grid.isIndexValid(index))
            {
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double coeff = LaplacainCoefficients(ManhanDistance,h);
                tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
            }
            //ghost cell
            else
            {
                //get the coffe of DiffOperator
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double diff_coeff = LaplacainCoefficients(ManhanDistance,h);
                //get the coeff of GhostCell
                int layer = grid.layer(index);
                Normal inner_normal = centered_cell - index;
                Normal direction = normalize(inner_normal);
                //start from the boundary cell
                MultiIndex boundary_cell = index + direction*layer;
                for(int k=0; k<3; k++)
                {
                    double ghost_coeff = GhostCellCoefficients(layer, k, h, boundary_type);
                    tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(boundary_cell+direction*k), ghost_coeff*diff_coeff));
                }
                rhs(centered_cell) += diff_coeff*GhostCellCoefficients(layer, -1, h, boundary_type)*boudary_value(boundary_cell, -1*direction);
            }
        }
    }
    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}


/** 
DivergenceOperator::DivergenceOperator(const Grid &grid, const VectorBoundaryFaceAvr &boudary_value, type boundary_type): rhs(grid) 
{
    int n_cells = grid.get_size()[0] * grid.get_size()[1];
    double h = grid.get_h();
    TripletList tripletList;
    matrix.resize(n_cells, n_cells);

    loop_inner_cell_2(grid,i,j)
    {
        MultiIndex centered_cell{i,j};
        MultiIndexList index_diff_cells_x = FourthOrderDiffCells_X(centered_cell);
        MultiIndexList index_diff_cells_y = FourthOrderDiffCells_Y(centered_cell);
        for (const auto &index : index_diff_cells_x)
        {
            int ManhanDistance = SignManhattanDistance(centered_cell, index);
            double coeff = GradientCoefficients(ManhanDistance,h);
            tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
        }
        for (const auto &index : index_diff_cells_y)
        {
            int ManhanDistance = SignManhattanDistance(centered_cell, index);
            double coeff = GradientCoefficients(ManhanDistance,h);
            tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
        }
    }
    MultiIndexSet outer_cells = getOuterCells(grid);
    for (const auto &centered_cell : outer_cells)
    {
        MultiIndexList index_diff_cells_x = FourthOrderDiffCells_X(centered_cell);
        MultiIndexList index_diff_cells_y = FourthOrderDiffCells_Y(centered_cell);
        for (const auto &index : index_diff_cells_x)
        {
            if (grid.isIndexValid(index))
            {
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double coeff = GradientCoefficients(ManhanDistance,h);
                tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
            }
            //ghost cell
            else
            {
                //get the coffe of DiffOperator
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double diff_coeff = GradientCoefficients(ManhanDistance,h);
                //get the coeff of GhostCell
                int layer = grid.layer(index);
                Normal inner_normal = centered_cell - index;
                Normal direction = normalize(inner_normal);
                //start from the boundary cell
                MultiIndex boundary_cell = index + direction*layer;
                for(int k=0; k<3; k++)
                {
                    double ghost_coeff = GhostCellCoefficients(layer, k, h, boundary_type);
                    tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(boundary_cell+direction*k), ghost_coeff*diff_coeff));
                }
                rhs(centered_cell) += diff_coeff*GhostCellCoefficients(layer, -1, h, boundary_type)*boudary_value(boundary_cell, -1*direction)[0];
            }
        }
        for (const auto &index : index_diff_cells_y)
        {
            if (grid.isIndexValid(index))
            {
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double coeff = GradientCoefficients(ManhanDistance,h);
                tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(index), coeff));
            }
            //ghost cell
            else
            {
                //get the coffe of DiffOperator
                int ManhanDistance = SignManhattanDistance(centered_cell, index);
                double diff_coeff = GradientCoefficients(ManhanDistance,h);
                //get the coeff of GhostCell
                int layer = grid.layer(index);
                Normal inner_normal = centered_cell - index;
                Normal direction = normalize(inner_normal);
                //start from the boundary cell
                MultiIndex boundary_cell = index + direction*layer;
                for(int k=0; k<3; k++)
                {
                    double ghost_coeff = GhostCellCoefficients(layer, k, h, boundary_type);
                    tripletList.push_back(Triplet(grid.MultiToSingle(centered_cell), grid.MultiToSingle(boundary_cell+direction*k), ghost_coeff*diff_coeff));
                }
                rhs(centered_cell) += diff_coeff*GhostCellCoefficients(layer, -1, h, boundary_type)*boudary_value(boundary_cell, -1*direction)[1];
            }
        }
    }
    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}
*/