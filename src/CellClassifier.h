#pragma once

#include <iostream>
#include "../src/Vec.h"
#include "../src/Grid.h"
#include "../src/BoundaryCurve.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>

/// @brief Hash function for MultiIndex to be used in unordered_set 
namespace std {
    template <>
    struct hash<MultiIndex> {
        std::size_t operator()(const MultiIndex& k) const {
            return std::hash<int>()(k[0]) ^ (std::hash<int>()(k[1]) << 1);
        }
    };
}

using MultiIndexSet = std::unordered_set<MultiIndex>;
using VecList = std::vector<Vec>;
using CutPointsMap = std::unordered_map<MultiIndex, VecList>;

/// @brief Classify the cells in the grid into different categories
class CellClassifier
{
public:
    CellClassifier() = default;
    CellClassifier(const Grid &grid, const BoundaryCurve &boundaryCurve);
    
    const MultiIndexSet& getDeadCells() const { return DeadCells; }
    const MultiIndexSet& getAliveCells() const { return AliveCells; }
    const MultiIndexSet& getCutCells() const { return CutCells; }
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
    CutPointsMap CutPoints;

    MultiIndexSet SideCells;
    MultiIndexSet EdgeCells;
    MultiIndexSet CoreCells;
};

void CellClassifier::LocateCutCells(const Grid &grid, const BoundaryCurve &boundaryCurve)
{
    for(const auto& point : boundaryCurve.GetControlPoints())
    {
        MultiIndex index = grid.LocateCell(point);
        CutCells.insert(index);
        CutPoints[index].push_back(point);
    }
}
