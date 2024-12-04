#pragma once

#include <iostream>
#include "../src/Vec.h"
#include "../src/Grid.h"
#include "../src/BoundaryCurve.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <iomanip>

/// @brief Hash function for MultiIndex to be used in unordered_set 
namespace std {
    template <>
    struct hash<MultiIndex> {
        std::size_t operator()(const MultiIndex& k) const {
            return std::hash<int>()(k[0]) ^ (std::hash<int>()(k[1]) << 1);
        }
    };
}
/// @brief Rename
using ParaInterval = std::pair<double, double>;
using MultiIndex = std::array<int, 2>;

using ParaIntervalList = std::vector<ParaInterval>;
using MultiIndexList = std::vector<MultiIndex>;
using MultiIndexSet = std::unordered_set<MultiIndex>;


/// @brief Given the spline, and the direction of curve passing through the cell
/// @param grid 
/// @param spline spline
/// @param interval <a,b> the interval of the curve in the cell
/// @param normal outer normal of the face
/// @param physical_tol physical distance tolerance
/// @param para_tol parameter interval length tolerance 
/// @param max_iter maximum iteration times
/// @return bisection parameter
double IntervalBisection(const Grid &grid, const Spline2d &spline, const ParaInterval &interval, MultiIndex cellindex, Normal normal, double physical_tol, double para_tol, int max_iter);
/// @brief Compute the difference of two MultiIndex to get the normal of the face
/// @param a 
/// @param b 
/// @return outer normal of the face
Normal minus(const MultiIndex& a, const MultiIndex& b);

/// @brief Compute the which piece of curve belonging to which cell
class CurveBelonging
{
public:
    CurveBelonging() = default;


    /// @brief Roughly check which cell is cut by the curve
    void RoughlyCheck(const Grid &grid, const Spline2d &spline, double step = 0.1);
    /// @brief Check if the curve spans two cells
    bool UncontinuousCell();
    /// @brief Adaptive check which cell is cut by the curve
    void AdaptiveCheck(const Grid &grid, const Spline2d &spline);
    /// @brief Computering the intersection of the curve and the cell
    //Currently using the bisection, Newton's method will be added in the future
    void PieceWiseBelonging(const Grid &grid, const Spline2d &spline, double physical_tol, double para_tol, int max_iter);
    VecList getIntersectionPoints(const Spline2d &spline) const;
    
    

    const ParaIntervalList& getParaIntervals() const { return ParaIntervals; }
    const MultiIndexList& getMultiIndices() const { return MultiIndices; }
    const MultiIndexSet& getLocalCutCells() const { return LocalCutCells; }

protected:
    ParaIntervalList ParaIntervals;
    MultiIndexList MultiIndices;
    MultiIndexSet LocalCutCells;
};
/******Function Implementation*****/
/**********************************/
void CurveBelonging::RoughlyCheck(const Grid &grid, const Spline2d &spline, double step)
{
    ParaIntervals.clear();
    MultiIndices.clear();
    double start = 0.0;
    double end = 1.0;

    MultiIndex index_first = grid.LocateCell(Vec{spline(start)});
    LocalCutCells.insert(index_first);
    MultiIndices.push_back(index_first);
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
        double last = ParaIntervals.back().second;
        if (last < 1.0)
        {
            ParaIntervals.push_back({last, 1.0});
        }
    }
    else
    {
        ParaIntervals.push_back({0.0, 1.0});
    }
}
bool CurveBelonging::UncontinuousCell()
{
    bool uncontinuous = false;
    for(size_t i = 1; i < MultiIndices.size(); ++i)
    {
        int diff_i = std::abs(MultiIndices[i][0] - MultiIndices[i - 1][0]);
        int diff_j = std::abs(MultiIndices[i][1] - MultiIndices[i - 1][1]);
        int diff = diff_i + diff_j;
        if(diff > 1)
        {
            uncontinuous = true;
            break;
        }
    }
    return uncontinuous;
}
void CurveBelonging::AdaptiveCheck(const Grid &grid, const Spline2d &spline)
{
    double step = 0.1;
    RoughlyCheck(grid, spline, step);
    while(UncontinuousCell())
    {
        step /= 2;
        RoughlyCheck(grid, spline, step);
    }
}

void CurveBelonging::PieceWiseBelonging(const Grid &grid, const Spline2d &spline, double physical_tol, double para_tol, int max_iter)
{
    int n_cutcell = LocalCutCells.size();
    for(int i=0; i<n_cutcell-1; i++)
    {
        MultiIndex current_cell = MultiIndices[i];
        MultiIndex next_cell = MultiIndices[i+1];
        Normal normal = minus(next_cell, current_cell);
        ParaInterval interval = ParaIntervals[i];
        double bisection_result = IntervalBisection(grid, spline, interval, current_cell, normal, physical_tol, para_tol, max_iter);
        ParaIntervals[i].second = bisection_result;
    }
    for(int i=1; i<n_cutcell; i++)
    {
        ParaIntervals[i].first = ParaIntervals[i-1].second;
    }
}
VecList CurveBelonging::getIntersectionPoints(const Spline2d &spline) const
{
    VecList intersection_list;
    for(const auto& interval : ParaIntervals)
    {
        intersection_list.push_back(Vec{spline(interval.second)});
    }
    intersection_list.pop_back();
    return intersection_list;
}

/*******************************************************************************************/
double IntervalBisection(const Grid &grid, const Spline2d &spline, const ParaInterval &interval, MultiIndex cellindex, Normal normal, double physical_tol, double para_tol, int max_iter = 100)
{
    //parameters of curve
    double a = interval.first;
    double b = interval.second;
    double mid = (a + b) / 2;
    //end points of the interval
    Vec left_end{ spline(a) };
    Vec right_end{ spline(b) };
    //mid point of the interval
    Vec mid_point{ spline(mid) };
    //sign distance of the end points and mid point to the face
    double mid_dis = grid.SignDistance(cellindex, normal, mid_point);
    double left_dis = grid.SignDistance(cellindex, normal, left_end);
    double right_dis = grid.SignDistance(cellindex, normal, right_end);
    if(left_dis * right_dis > 0)
    {
        std::cerr << "The interval is not cut by the curve" << std::endl;
        return -1;
    }
    if(left_dis==0)
    {
        return a;
    }
    if (right_dis == 0)
    {
        return b;
    }
    while (std::fabs(mid_dis) >= physical_tol && (b - a) > para_tol && max_iter-- > 0)
    {   
        assert(left_dis * right_dis <= 0 || right_dis * mid_dis <= 0);

        if (left_dis * mid_dis < 0)
        {
            b = mid;
        }
        else
        {
            a = mid;
        }

        mid = (a + b) / 2;

        mid_point = Vec{spline(mid)};
        left_end = Vec{spline(a)};
        right_end = Vec{spline(b)};

        mid_dis = grid.SignDistance(cellindex, normal, mid_point);
        left_dis = grid.SignDistance(cellindex, normal, left_end);
        right_dis = grid.SignDistance(cellindex, normal, right_end);
    }
    return mid;
}

Normal minus(const MultiIndex& a, const MultiIndex& b)
{
    return {a[0] - b[0], a[1] - b[1]};
}