#pragma once

#include <iostream>
#include "../../lib/Eigen/Core"
#include "../../lib/unsupported/Eigen/Splines"
#include "../Vec.h"
#include "../Grid.h"
#include "../MultiIndexSet.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <iomanip>

/// @brief Rename
using Spline2d = Eigen::Spline<double, 2, 3>;
using VecList = std::vector<Vec>;

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





/// @brief Compute the which piece of curve belonging to which cell
class CurveBelonging
{
public:
    CurveBelonging() = default;


    /// @brief Roughly check which cell is cut by the curve
    void RoughlyCheck(const Grid &grid, const Spline2d &spline, double step);
    /// @brief Check if the curve spans two cells
    bool UncontinuousCell();
    /// @brief Adaptive check which cell is cut by the curve
    void AdaptiveCheck(const Grid &grid, const Spline2d &spline);
    /// @brief Computering the intersection of the curve and the cell
    //Currently using the bisection, Newton's method will be added in the future
    void PieceWiseBelonging(const Grid &grid, const Spline2d &spline, double physical_tol, double para_tol, int max_iter);
    /// @brief Get the intersection points of the curve and the cell. 
    /// If the curve locates in a single cell, the function will return  A EMPTY VECTOR
    /// @param spline
    /// @return std::vector<Vec> intersection points
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
    //Locating the first cell
    LocalCutCells.insert(index_first);
    MultiIndices.push_back(index_first);

    //Roughly check which cell is cut by the curve
    //If a cell is cut by the curve, then the cell is added to the MultiIndices
    //And correspanding interval will be stored
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
    //The last interval should be added, the end piece whole lies in the last cell
    if (!ParaIntervals.empty())
    {
        double last = ParaIntervals.back().second;
        if (last < end)
        {
            ParaIntervals.push_back({last, end});
        }
    }
    //If there is not any cell cut by the curve, then the curve is in a single cell
    else
    {
        ParaIntervals.push_back({start, end});
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
        //If the curve spans two cells, we name it uncontinuous
        //The adjacent cells that are stored should be adjacent in physical space.
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
    int n_cutcell = MultiIndices.size();
    //As shown before, the last piece of curve lies in a single cell
    for(size_t i = 0; i < n_cutcell-1; i++)
    {
        MultiIndex current_cell = MultiIndices[i];
        MultiIndex next_cell = MultiIndices[i+1];
        //The adjacent cells tell the direction of the curve
        Normal normal = next_cell - current_cell;
        ParaInterval interval = ParaIntervals[i];
        double bisection_result = IntervalBisection(grid, spline, interval, current_cell, normal, physical_tol, para_tol, max_iter);
        ParaIntervals[i].second = bisection_result;
    }
    for(size_t i = 1; i < n_cutcell; i++)
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

// Normal minus(const MultiIndex& a, const MultiIndex& b)
// {
//     return {a[0] - b[0], a[1] - b[1]};
// }

