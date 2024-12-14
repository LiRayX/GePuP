#pragma once


#include "Vec.h"
#include "Grid.h"
#include "CyclicCurve.h"
#include "MultiIndexSet.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <iostream>

using VecList = std::vector<Vec>;
using CutCellMapping = std::unordered_map<MultiIndex, std::set<double>>;


class CellDivision
{
public:
    void LocateDeadCells(const Grid &grid);
    void LocateAllCells(const Grid &grid, const CyclicCurve &boundaryCurve);
    void LocateCutCells(const Grid &grid, const CyclicCurve &boundaryCurve);


    void LocateSideCells(const Grid &grid);
    // void LocateEdgeCells(const Grid &grid);

    void LocateCoreCells(const Grid &grid);
    bool isCoreCell(const MultiIndex &index, const Grid &grid);
        
    void LocateExtendedCells(const Grid &grid);
    bool isExtendedCell(const MultiIndex &index, const Grid &grid);

    const MultiIndexSet& getDeadCells() const { return DeadCells; }
    const MultiIndexSet& getAliveCells() const { return AliveCells; }
    const MultiIndexSet& getCutCells() const { return CutCells; }
    const MultiIndexSet& getSideCells() const { return SideCells; }
    const MultiIndexSet& getEdgeCells() const { return EdgeCells; }
    const MultiIndexSet& getCoreCells() const { return CoreCells; }
    const MultiIndexSet& getExtendedCells() const { return ExtendedCells; }


    const CutCellMapping& getCutCellInfo() const { return CutCellInfo; }
protected:

    MultiIndexSet DeadCells;
    MultiIndexSet AliveCells;

    MultiIndexSet CutCells;
    CutCellMapping CutCellInfo;

    MultiIndexSet SideCells;

    MultiIndexSet EdgeCells;
    MultiIndexSet CoreCells;

    MultiIndexSet ExtendedCells;
};
void CellDivision::LocateAllCells(const Grid &grid, const CyclicCurve &boundaryCurve)
{
    this->LocateCutCells(grid, boundaryCurve);
    this->LocateDeadCells(grid);
    this->LocateSideCells(grid);
    // this->LocateEdgeCells(grid);
    this->LocateCoreCells(grid);
    this->LocateExtendedCells(grid);
}



void CellDivision::LocateCutCells(const Grid &grid, const CyclicCurve &boundaryCycle)
{
  ParaSet intersections = boundaryCycle.getIntersections();
  for (const auto& para: intersections)
  {
    Vec point = boundaryCycle.getPoint(para);
    MultiIndex index = grid.LocateCell(point);
    Normal normal = grid.getDirection(index, point, 1e-6);
    CutCells.insert(index);
    CutCells.insert(index + normal);
    CutCellInfo[index].insert(para);
    CutCellInfo[index + normal].insert(para);
  }
}

void CellDivision::LocateDeadCells(const Grid &grid)
{
    RangeMap rangeMap = getRangeOfCutCells(CutCells);
    for(const auto &map : rangeMap)
    {
        int i = map.first;
        int j_min = map.second.first;
        int j_max = map.second.second;
        for(int j = j_min+1; j <= j_max-1; j++)
        {
            if(CutCells.find({i, j}) == CutCells.end())
            {
                DeadCells.insert({i, j});
            }
        }
    }
}
void CellDivision::LocateSideCells(const Grid &grid)
{
    for (const auto& index : CutCells)
    {
        std::vector<MultiIndex> neighbours = VonNeumannNeighbour(index, grid);
        for (const auto& neighbour : neighbours)
        {
            //If not a cut cell or a dead cell
            if (CutCells.find(neighbour) == CutCells.end() && DeadCells.find(neighbour) == DeadCells.end())
            {
                SideCells.insert(neighbour);
            }
        }
    }
}

/// @brief Locate the core cells and edge cells
void CellDivision::LocateCoreCells(const Grid &grid)
{
    //Dead cells may contain in the 2nd-inner part of the domain
    loop_inner_cell_2(grid, i, j)
    {
        MultiIndex index{i, j};
        //Excluding cut cells, dead cells, side cells.
        if (CutCells.find(index) == CutCells.end() && DeadCells.find(index) == DeadCells.end() && SideCells.find(index) == SideCells.end())
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
void CellDivision::LocateExtendedCells(const Grid &grid)
{
    MultiIndexSet OuterCells = getOuterCells(grid);
    for (const auto& index : OuterCells)
    {
        //Excluding cut cells, dead cells, side cells.
        if (CutCells.find(index) == CutCells.end() && DeadCells.find(index) == DeadCells.end() && SideCells.find(index) == SideCells.end())
        {
            if (isExtendedCell(index, grid))
            {
                ExtendedCells.insert(index);
            }
            else
            {
                EdgeCells.insert(index);
            }
        }
    }
}


/// @brief True means core cell, False means Edge cell
/// @param index 
/// @param grid 
/// @return bool
bool CellDivision::isCoreCell(const MultiIndex &index, const Grid &grid)
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
/// @brief True means core cell, False means Edge cell
/// @param index 
/// @param grid 
/// @return bool
bool CellDivision::isExtendedCell(const MultiIndex &index, const Grid &grid)
{
    bool isExtend = true;
    std::vector<MultiIndex> neighbours = ExtendedVonNeumannNeighbour(index, grid);
    for (const auto& neighbour : neighbours)
    {
        if (CutCells.find(neighbour) != CutCells.end() || DeadCells.find(neighbour) != DeadCells.end())
        {
            isExtend = false;
            break;
        }
    }
    return isExtend;
}

// /******************************************************************************* */
// /// @brief Dropped
// /// @param grid 
// void CellDivision::LocateEdgeCells(const Grid &grid)
// {
//     for (const auto& index : SideCells)
//     {
//         std::vector<MultiIndex> neighbours = VonNeumannNeighbour(index, grid);
//         for (const auto& neighbour : neighbours)
//         {
//             //If not a cut cell or a dead cell or a side cell
//             if (CutCells.find(neighbour) == CutCells.end() && DeadCells.find(neighbour) == DeadCells.end() && SideCells.find(neighbour) == SideCells.end())
//             {
//                 EdgeCells.insert(neighbour);
//             }
//         }
//     }
// }