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

using ParaInterval = std::pair<double, double>;
using MultiIndex = std::array<int, 2>;

using ParaIntervalList = std::vector<ParaInterval>;
using MultiIndexList = std::vector<MultiIndex>;
using MultiIndexSet = std::unordered_set<MultiIndex>;


/// @brief Compute the which piece of curve belonging to which cell
class CurveBelonging
{
public:
    CurveBelonging() = default;

    void Bisection(const Grid &grid, const Spline2d &spline, double step = 0.1, double tol = 1e-6);

    const ParaIntervalList& getParaIntervals() const { return ParaIntervals; }
    const MultiIndexList& getMultiIndices() const { return MultiIndices; }
    const MultiIndexSet& getLocalCutCells() const { return LocalCutCells; }

protected:
    ParaIntervalList ParaIntervals;
    MultiIndexList MultiIndices;
    MultiIndexSet LocalCutCells;
};

void CurveBelonging::Bisection(const Grid &grid, const Spline2d &spline, double step, double tol)
{
    double h = grid.get_h();

    double start = 0.0;
    double end = 1.0;

    MultiIndex index_first = grid.LocateCell(Vec{spline(0)});

    //Roughly check which cell is cut by the curve
    for (double u = step; u <= end; u += step)
    {
        Vec point{spline(u)};
        MultiIndex index_second = grid.LocateCell(point);
        if (index_second != index_first)
        {
            LocalCutCells.insert(index_second);
            ParaIntervals.push_back({start, u});
            MultiIndices.push_back(index_second);
            start = u;
            index_first = index_second;
        }
    }
    if (!ParaIntervals.empty())
    {
        ParaIntervals.back().second = 1.0;
    }
    else
    {
        ParaIntervals.push_back({0.0, 1.0});
        MultiIndices.push_back(index_first);
    }

    // MultiIndex index_start = grid.LocateCell(Vec{spline(0)});
    // MultiIndex index_end = grid.LocateCell(Vec{spline(1)});
    // int n_cutcell = std::abs(index_end[0] - index_start[0]) + std::abs(index_end[1] - index_start[1]);
    // if (index_start == index_end)
    // {
    //     ParaIntervals.push_back({start, end});
    //     MultiIndices.push_back(index_start);
    // } 
    // else
    // { 
    //     while(in)

    // }
}
