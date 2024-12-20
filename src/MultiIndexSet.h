#pragma once

#include <unordered_set>
#include <unordered_map>
#include <limits>
#include "../lib/Eigen/Dense"
#include "Grid.h"

using Normal = std::array<int, 2>;
using MultiIndex = std::array<int, 2>;
using MultiIndexSet = std::unordered_set<MultiIndex>;
using RangeMap = std::unordered_map<int, std::pair<int, int>>;


/// @brief Compute the difference of two MultiIndex to get the normal of the face
Normal operator-(const MultiIndex& a, const MultiIndex& b);
/// @brief Find next cell index
MultiIndex operator+(const MultiIndex& a, const Normal& b);
/// @brief VonNeumann Neighbour
std::vector<MultiIndex> VonNeumannNeighbour(const MultiIndex &index);
std::vector<MultiIndex> VonNeumannNeighbour(const MultiIndex &index, const Grid &grid);
std::vector<MultiIndex> ExtendedVonNeumannNeighbour(const MultiIndex &index, const Grid &grid);
/// @brief Neighbors for Ghost Cell
std::vector<MultiIndex> GhostNeighbour(const MultiIndex &index, const Grid &grid);
/// @brief Binary Operator of sets.
MultiIndexSet unionSets(const MultiIndexSet& set1, const MultiIndexSet& set2);
MultiIndexSet differenceSets(const MultiIndexSet& set1, const MultiIndexSet& set2); 

MultiIndexSet getOuterCells(const Grid& grid);

std::ostream &operator<<(std::ostream &os, const MultiIndex &mi);

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


Normal operator-(const MultiIndex& a, const MultiIndex& b)
{
    return {a[0] - b[0], a[1] - b[1]};
}

MultiIndex operator+(const MultiIndex& a, const Normal& b)
{
    return {a[0] + b[0], a[1] + b[1]};
}
std::vector<MultiIndex> VonNeumannNeighbour(const MultiIndex &index)
{
    std::vector<MultiIndex> neighbour; 
    std::vector<Normal> normals = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};  
    for (const auto& normal : normals) 
    {
        neighbour.push_back(index + normal);
    }
    return neighbour;
}

std::vector<MultiIndex> VonNeumannNeighbour(const MultiIndex &index, const Grid &grid)
{
    std::vector<MultiIndex> neighbour; 
    std::vector<Normal> normals = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

    for (const auto& normal : normals) 
    {
        MultiIndex neighbourIndex = index + normal;
        if (grid.isIndexValid(neighbourIndex)) 
        {
            neighbour.push_back(neighbourIndex);
        }
    }
    return neighbour;
}
std::vector<MultiIndex> ExtendedVonNeumannNeighbour(const MultiIndex &index, const Grid &grid)
{
    std::vector<MultiIndex> neighbour; 
    std::vector<Normal> normals = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}, {0, 2}, {0, -2}, {2, 0}, {-2, 0}};
    for (const auto& normal : normals) 
    {
        MultiIndex neighbourIndex = index + normal;
        if (grid.isIndexValid(neighbourIndex)) 
        {
            neighbour.push_back(neighbourIndex);
        }
    }
    return neighbour;
}


std::vector<MultiIndex> GhostNeighbour(const MultiIndex &index, const Grid &grid)
{
    std::vector<MultiIndex> neighbour; 

    int i = index[0];
    int j = index[1];

    int m = grid.get_size()[0];
    int n = grid.get_size()[1];

    int n_ghost = 4;
    //Left
    if(i == 0 || i ==1)
    {
        for(int k = 0; k < n_ghost; k++)
        {
            neighbour.push_back({k, j});
        }
    }
    //Right
    if(i == m-1 || i == m-2)
    {
        for(int k = m-1; k > m-1-n_ghost; k--)
        {
            neighbour.push_back({k, j});
        }
    }
    //Down
    if(j == 0 || j == 1)
    {
        for(int k = 0; k < n_ghost; k++)
        {
            neighbour.push_back({i, k});
        }
    }
    //Up
    if(j == n-1 || j == n-2)
    {
        for(int k = n-1; k > n-1-n_ghost; k--)
        {
            neighbour.push_back({i, k});
        }
    }
    return neighbour;
}



// Union of two sets
MultiIndexSet unionSets(const MultiIndexSet& set1, const MultiIndexSet& set2) 
{
    MultiIndexSet result = set1;
    result.insert(set2.begin(), set2.end());
    return result;
}

// Differ of two sets
MultiIndexSet differenceSets(const MultiIndexSet& set1, const MultiIndexSet& set2) 
{
    MultiIndexSet result;
    for (const auto& mi : set1) {
        if (set2.find(mi) == set2.end()) {
            result.insert(mi);
        }
    }
    return result;
}
/// @brief Get the outer cells of the grid, which are the cells on the boundary or near the boundary need a ghost cell.
/// @param grid 
/// @return 
MultiIndexSet getOuterCells(const Grid& gd)
{
    MultiIndexSet outerCells;
    for (int i1 = 0; i1 < 2; i1++)                  
        for (int i0 = 0; i0 < gd.get_size()[0]; i0++)
            outerCells.insert({i0, i1}); 
    for (int i1 = gd.get_size()[1]-2; i1 < gd.get_size()[1]; i1++) 
        for (int i0 = 0; i0 < gd.get_size()[0]; i0++) 
            outerCells.insert({i0, i1});
    for (int i1 = 0; i1 < 2; i1++)                  
        for (int i0 = 2; i0 < gd.get_size()[1]-3; i1++)
            outerCells.insert({i0, i1});
    for (int i1 = gd.get_size()[1]-2; i1 < gd.get_size()[1]; i1++) 
        for (int i0 = 2; i0 < gd.get_size()[0]-2; i0++)
            outerCells.insert({i0, i1});
    // int m = grid.get_size()[0];
    // int n = grid.get_size()[1];
    // loop_cell_2(grid, i, j)
    // {
    //     if(i == 0 || i == m-1 || j == 0 || j == n-1 || i == 1 || i == m-2 || j == 1 || j == n-2)
    //     {
    //         outerCells.insert({i, j});
    //     }
    // }
    return outerCells;
}

std::ostream &operator<<(std::ostream &os, const MultiIndex &mi)
{
    os << "(" << mi[0] << ", " << mi[1] << ")";
    return os;
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


