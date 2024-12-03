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



    void RoughlyCheck(const Grid &grid, const Spline2d &spline, double step = 0.1);
    void AdaptiveCheck(const Grid &grid, const Spline2d &spline);
    void Bisection(const Grid &grid, const Spline2d &spline, double tol = 1e-6);

    const ParaIntervalList& getParaIntervals() const { return ParaIntervals; }
    const MultiIndexList& getMultiIndices() const { return MultiIndices; }
    const MultiIndexSet& getLocalCutCells() const { return LocalCutCells; }

protected:
    ParaIntervalList ParaIntervals;
    MultiIndexList MultiIndices;
    MultiIndexSet LocalCutCells;
};
void CurveBelonging::RoughlyCheck(const Grid &grid, const Spline2d &spline, double step)
{
    ParaIntervals.clear();
    MultiIndices.clear();
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
}

void CurveBelonging::AdaptiveCheck(const Grid &grid, const Spline2d &spline)
{
    double step = 0.1;
    RoughlyCheck(grid, spline, step);
     for (size_t k = 1; k < MultiIndices.size(); ++k) 
     {
        int diff_i = std::abs(MultiIndices[k][0] - MultiIndices[k - 1][0]);
        int diff_j = std::abs(MultiIndices[k][1] - MultiIndices[k - 1][1]);

        // 如果 i 或 j 指标之差超过了 1，则对对应位置的参数进行加细
        if (diff_i + diff_j > 1) 
        {
            RoughlyCheck(grid, spline, step / 2);
        }
    }

}

void CurveBelonging::Bisection(const Grid &grid, const Spline2d &spline, double tol)
{
    
}
