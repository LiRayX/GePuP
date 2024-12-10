#pragma once

#include <unordered_set>
#include <unordered_map>
#include <limits>
#include "../lib/Eigen/Dense"
#include "Grid.h"

using MultiIndex = std::array<int, 2>;
using MultiIndexSet = std::unordered_set<MultiIndex>;
using RangeMap = std::unordered_map<int, std::pair<int, int>>;

/// @brief Hash function for MultiIndex to be used in unordered_set 
namespace std 
{
    template <>
    struct hash<MultiIndex> 
    {
        std::size_t operator()(const MultiIndex& k) const 
        {
            return std::hash<int>()(k[0]) ^ (std::hash<int>()(k[1]) << 1);
        }
    };
}
/// @brief Get the range of the first index of the MultiIndexSet
/// @param CutCells 
/// @param grid 
/// @return 
RangeMap getRangeOfCutCells(const MultiIndexSet& CutCells) 
{
    RangeMap rangeMap;
    for (const auto& cell : CutCells) {
        int i = cell[0];
        int j = cell[1];
        if (rangeMap.find(i) == rangeMap.end()) 
        {
            rangeMap[i] = {j, j};
        } 
        else 
        {
            if (j < rangeMap[i].first) 
            {
                rangeMap[i].first = j;
            }
            if (j > rangeMap[i].second) 
            {
                rangeMap[i].second = j;
            }
        }
    }
    return rangeMap;
}

/// @brief Matrix to store the cut cells
/// @param multiIndexSet 
/// @param grid 
/// @return 
Eigen::MatrixXi convertToEigenMatrix(const MultiIndexSet& multiIndexSet, const Grid& grid)
{
    int rows = grid.get_size()[0];
    int cols = grid.get_size()[1];
    Eigen::MatrixXi matrix = Eigen::MatrixXi::Constant(rows, cols, -1);
    //Set cutcellls to 0.
    for (const auto& index : multiIndexSet) {
        int i = index[0];
        int j = index[1];
        if (i >= 0 && i < rows && j >= 0 && j < cols) 
        {
            matrix(i, j) = 0;
        }
    }
    return matrix;
}




/// @brief Get the range of the first index of the MultiIndexSet
/// @param CutCells 
/// @return 
std::pair<int, int> getRangeOfFirstIndex(const MultiIndexSet& CutCells) 
{
    int minIndex = std::numeric_limits<int>::max();
    int maxIndex = std::numeric_limits<int>::min();

    for (const auto& cell : CutCells) {
        if (cell[0] < minIndex) {
            minIndex = cell[0];
        }
        if (cell[0] > maxIndex) {
            maxIndex = cell[0];
        }
    }
    return {minIndex, maxIndex};
}

/// @brief Given a first index, get the range of the second index
std::pair<int, int> getRangeOfSecondIndex(const MultiIndexSet& CutCells, int firstIndex) 
{
    int minIndex = std::numeric_limits<int>::max();
    int maxIndex = std::numeric_limits<int>::min();

    for (const auto& cell : CutCells) {
        if (cell[0] == firstIndex) {
            if (cell[1] < minIndex) {
                minIndex = cell[1];
            }
            if (cell[1] > maxIndex) {
                maxIndex = cell[1];
            }
        }
    }
    return {minIndex, maxIndex};
}