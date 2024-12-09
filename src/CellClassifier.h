#pragma once

#include <iostream>
#include "Vec.h"
#include "Grid.h"
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

using MultiIndexSet = std::unordered_set<MultiIndex>;
using VecList = std::vector<Vec>;
using CutCellMap = std::unordered_map<MultiIndex, CutCellInfo>;


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

    void LocateCutCells(const Grid &grid, const BoundaryCurve &boundary);
    void LocateSideCells(const Grid &grid, const BoundaryCurve &boundary);
    void LocateEdgeCells(const Grid &grid, const BoundaryCurve &boundary);
    void LocateCoreCells(const Grid &grid, const BoundaryCurve &boundary);
protected:
    MultiIndexSet DeadCells;
    MultiIndexSet AliveCells;

    MultiIndexSet CutCells;
    CutCellMap CutCellInfo;

    MultiIndexSet SideCells;
    MultiIndexSet EdgeCells;
    MultiIndexSet CoreCells;
};

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



// void CellClassifier::LocateCutCells(const Grid &grid, const BoundaryCurve &boundaryCurve)
// {
//     for(const auto& point : boundaryCurve.GetControlPoints())
//     {
//         MultiIndex index = grid.LocateCell(point);
//         CutCells.insert(index);
//         CutPoints[index].push_back(point);
//     }
// }


