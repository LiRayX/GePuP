#pragma once

#include <iostream>
#include "../Vec.h"
#include "../Grid.h"
#include "BoundaryCurve.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>


struct CutCellInfo
{
    std::vector<int> index_spline;
    std::vector<ParaInterval> parainterval;
};


using VecList = std::vector<Vec>;
using CutCellMap = std::unordered_map<MultiIndex, CutCellInfo>;


// /// @brief Compute the normal_average of the curve in the cell
// Vec getNormal_CutCell(const MultiIndex &index, const CutCellInfo &info, const BoundaryCurve &boundarycurve, double step);
// Vec getNormal_Sum(const Spline2d &spline, const ParaInterval &interval, double step);

/// @brief Classify the cells in the grid into different categories
class CellClassifier
{
public:
    CellClassifier() = default;
    CellClassifier(const Grid &grid, const BoundaryCurve &boundaryCurve);
    
    const MultiIndexSet& getDeadCells() const { return DeadCells; }
    const MultiIndexSet& getAliveCells() const { return AliveCells; }
    const MultiIndexSet& getCutCells() const { return CutCells; }
    const CutCellMap& getCutCellInfo() const { return CutCellInfo; }
    const MultiIndexSet& getSideCells() const { return SideCells; }
    const MultiIndexSet& getEdgeCells() const { return EdgeCells; }
    const MultiIndexSet& getCoreCells() const { return CoreCells; }

    void LocateDeadCells(const Grid &grid);
    void LocateAliveCells(const Grid &grid);
    void LocateCutCells(const Grid &grid, const BoundaryCurve &boundary);

    void LocateSideCells(const Grid &grid);
    bool isSideCell(const MultiIndex &index, const Grid &grid);

    // void LocateEdgeCells(const Grid &grid);
    // bool isEdgeCell(const MultiIndex &index, const Grid &grid);

    /// @brief Locate the core cells and edge cells
    void LocateCoreCells(const Grid &grid);
    /// @brief Verify if the cell is a core cell or an edge cell
    /// @param index 
    /// @param grid 
    /// @return true means core cell, false means edge cell
    bool isCoreCell(const MultiIndex &index, const Grid &grid);
protected:
    MultiIndexSet DeadCells;
    MultiIndexSet AliveCells;

    MultiIndexSet CutCells;
    CutCellMap CutCellInfo;

    MultiIndexSet SideCells;
    MultiIndexSet EdgeCells;
    MultiIndexSet CoreCells;
};

// /// @brief This function need LocateCutCells to be called first
// /// @param grid 
// /// @param boundaryCurve 
// void CellClassifier::LocateAliveCells(const Grid &grid)
// {
//     RangeMap rangeMap = getRangeOfCutCells(CutCells);
//     Eigen::MatrixXi indexMatrix = convertToEigenMatrix(CutCells, grid);
//     for(const auto &map : rangeMap)
//     {
//         int i = map.first;
//         int j_min = map.second.first;
//         int j_max = map.second.second;
//         for(int j = j_min; j <= j_max; j++)
//         {
//             AliveCells.insert({i, j});
//             if(indexMatrix(i, j-1) == 0 || indexMatrix(i, j+1) == 0 || indexMatrix(i-1, j) == 0 || indexMatrix(i+1, j) == 0)
//             {
//                 SideCells.insert({i, j});
//             }
//         }
//     }
// }

void CellClassifier::LocateAliveCells(const Grid &grid)
{
    this->LocateSideCells(grid);
    // this->LocateEdgeCells(grid);
    this->LocateCoreCells(grid);
}


void CellClassifier::LocateCutCells(const Grid &grid, const BoundaryCurve &boundaryCurve)
{
    const SplineList &splines = boundaryCurve.getSplines();
    int n_segement = splines.size();
    const PieceWiseBelongingList &PieceWiseBelongingIfo = boundaryCurve.getPieceWiseBelongingIfo();

    for(size_t i = 0; i < n_segement; i++)
    {
        //Get the i-th spline and its belonging information
        const Spline2d &spline = splines[i];
        const CurveBelonging &curveBelonging = PieceWiseBelongingIfo[i];
        
        //Cut cell index and its belonging parameters interval
        const MultiIndexList &multiIndices = curveBelonging.getMultiIndices();
        const ParaIntervalList &paraIntervals = curveBelonging.getParaIntervals();
        //Locate the cut cells and store the information
        for(size_t j = 0; j < multiIndices.size(); j++)
        {
            MultiIndex index = multiIndices[j];
            CutCells.insert(index);
            CutCellInfo[index].index_spline.push_back(i);
            CutCellInfo[index].parainterval.push_back(paraIntervals[j]);
        }
    }
}

void CellClassifier::LocateSideCells(const Grid &grid)
{
    RangeMap rangeMap = getRangeOfCutCells(CutCells);
    for(const auto &map : rangeMap)
    {
        int i = map.first;
        int j_min = map.second.first;
        int j_max = map.second.second;
        for(int j = j_min+1; j <= j_max-1; j++)
        {
            if(this->isSideCell({i, j}, grid))
            {
                SideCells.insert({i, j});
            }
        }
    }
}

void CellClassifier::LocateCoreCells(const Grid &grid)
{
    RangeMap rangeMap = getRangeOfCutCells(CutCells);
    for(const auto &map : rangeMap)
    {
        int i = map.first;
        int j_min = map.second.first;
        int j_max = map.second.second;
        for(int j = j_min+2; j <= j_max-2; j++)
        {
            MultiIndex index = {i, j};
            if(CutCells.find(index) != CutCells.end() || SideCells.find(index) != SideCells.end())
            {
                if (isCoreCell(index, grid))
                {
                    CoreCells.insert(index);
                }
                else
                {
                    EdgeCells.insert(index);
                }
            }
        }
    }
}


// void CellClassifier::LocateCoreCells(const Grid &grid)
// {
//     RangeMap rangeMap = getRangeOfCutCells(CutCells);
//     for(const auto &map : rangeMap)
//     {
//         int i = map.first;
//         int j_min = map.second.first;
//         int j_max = map.second.second;
//         for(int j = j_min+3; j <= j_max-3; j++)
//         {
//             if(this->isCoreCell({i, j}))
//             {
//                 CoreCells.insert({i, j});
//             }
//         }
//     }
// }







/// @brief Not cut cell, and at least one neighbour is a cut cell
/// @param index 
/// @return 
bool CellClassifier::isSideCell(const MultiIndex &index, const Grid &grid)
{
    //If the cell is a cut cell, then it is not a side cell
    if(CutCells.find(index) != CutCells.end())
    {
        return false;
    }
    std::vector<MultiIndex> neighbour = VonNeumannNeighbour(index, grid);
    for(const auto &neigh : neighbour)
    {
        if(CutCells.find(neigh) != CutCells.end())
        {
            return true;
        }
    }
    return false;
}



bool CellClassifier::isCoreCell(const MultiIndex &index, const Grid &grid)
{
    bool isCore = true;
    std::vector<MultiIndex> neighbours = ExtendedVonNeumannNeighbour(index, grid);
    for (const auto& neighbour : neighbours)
    {
        if (CutCells.find(neighbour) != CutCells.end() || DeadCells.find(neighbour) != DeadCells.end())
        {
            isCore = false;
            break;
        }
    }
    return isCore;
}

// /// @brief Not a cut cell or a side cell, and any neighbours is a side cell 
// /// @param index 
// /// @return 
// bool CellClassifier::isEdgeCell(const MultiIndex &index)
// {
//     //If the cell is a side cell, then it is not an edge cell
//     if(CutCells.find(index) != CutCells.end() || SideCells.find(index) != SideCells.end())
//     {
//         return false;
//     }
//     //If a neighbour of the cell is a side cell, and none of the neighbour is a cut cell, then the cell is an edge cell
//     std::vector<MultiIndex> neighbour = VonNeumannNeighbour(index);
//     for(const auto &neigh : neighbour)
//     {
//         if(SideCells.find(neigh) != SideCells.end())
//         {
//             return true;
//         }
//     }
//     return false;
// }

// bool CellClassifier::isCoreCell(const MultiIndex &index)
// {
//     if(CutCells.find(index) != CutCells.end() || SideCells.find(index) != SideCells.end() || EdgeCells.find(index) != EdgeCells.end())
//     {
//         return false;
//     }
//     else 
//     {
//         return true;
//     }
// }







// Vec getNormal_CutCell(const MultiIndex &index, const CutCellInfo &info, const BoundaryCurve &boundarycurve, double step)
// {
//     Vec normal{0.0, 0.0};

//     const std::vector<int> &index_spline = info.index_spline;
//     const std::vector<ParaInterval> &parainterval = info.parainterval;
//     const SplineList &splines = boundarycurve.getSplines();

//     int n_segment = index_spline.size();

//     for(int i = 0; i < n_segment; i++)
//     {
//         int index = index_spline[i];
//         const Spline2d &spline = splines[index];
//         const ParaInterval &interval = parainterval[i];
//         Vec current_normal = getNormal_Sum(spline, interval, step);
//         normal = normal + current_normal;
//     }
// }

// Vec getNormal_Sum(const Spline2d &spline, const ParaInterval &interval, double step)
// {
//     Vec normal{0.0, 0.0};
//     double start = interval.first;
//     double end = interval.second;
//     for(double u = start; u < end; u += step)
//     {
//        Vec current_normal = getNormal(spline, u);
//        normal = normal + current_normal;
//     }
//     return normal;
// }



// void CellClassifier::LocateCutCells(const Grid &grid, const BoundaryCurve &boundaryCurve)
// {
//     for(const auto& point : boundaryCurve.GetControlPoints())
//     {
//         MultiIndex index = grid.LocateCell(point);
//         CutCells.insert(index);
//         CutPoints[index].push_back(point);
//     }
// }


