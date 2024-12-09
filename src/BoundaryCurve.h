#pragma once


#include "../lib/Eigen/Core"
#include "../lib/unsupported/Eigen/Splines"
#include "Vec.h"
#include "Grid.h"
#include "CurveBelonging.h"
#include <iostream>
#include <vector>

using VecList = std::vector<Vec>;
using Spline2d = Eigen::Spline<double, 2, 3>;
using SplineList = std::vector<Spline2d>;
using PieceWiseBelongingList = std::vector<CurveBelonging>;


/// @brief Using segement B-spline to represent boundary curve
class BoundaryCurve
{
public:
    /// @brief Constructor
    BoundaryCurve() = default;
    BoundaryCurve(const VecList& controlPoints, const std::vector<int>& endPoints)
        : ControlPoints(controlPoints), IndexEndPoints(endPoints) {PieceWiseFit();}
    /// @brief Get ControlPoints
    VecList getControlPoints() const 
    { return ControlPoints; }
    /// @brief Get EndPoints
    const std::vector<int> &getEndPoints() const 
    { return IndexEndPoints; }
    /// @brief Get Segmented B-spline
    const SplineList &getSplines() const 
    { return Splines; }
    /// @brief Fit curve globally, a B-spline from startpoint to finalpoint
    /// @return Single B-Spline
    Spline2d GlobalFit() const;
    /// @brief Fit curve segemented, sequence B-splines according to IndexEndPoints
    void PieceWiseFit();
    
    void setPieceWiseBelongingIfo(const Grid &grid, double physical_tol, double para_tol, int max_iter);
    const PieceWiseBelongingList getPieceWiseBelongingIfo() const { return PieceWiseBelongingIfo; }


protected:
    VecList ControlPoints; //Start point should coincide with Final point, to get a closed curve
    std::vector<int> IndexEndPoints; //Deciding which points to be endpoints of piece-wise spline
    SplineList Splines; //Sequence of spline
    PieceWiseBelongingList PieceWiseBelongingIfo; //Asscoiated with each piece-wise spline
};



Spline2d BoundaryCurve::GlobalFit() const
{
    Eigen::MatrixXd pts(2, ControlPoints.size());
    for (int i = 0; i < ControlPoints.size(); ++i)
    {
        pts(0, i) = ControlPoints[i][0];
        pts(1, i) = ControlPoints[i][1];
    }
    return Eigen::SplineFitting<Spline2d>::Interpolate(pts, 3);
}

void BoundaryCurve::PieceWiseFit()
{
    for (int i = 0; i < IndexEndPoints.size() - 1; i++)
    {
        Eigen::MatrixXd pts(2, IndexEndPoints[i + 1] - IndexEndPoints[i] + 1);

        for (int j = IndexEndPoints[i], col = 0; j <= IndexEndPoints[i + 1]; ++j, ++col)
        {
            pts(0, col) = ControlPoints[j][0];
            pts(1, col) = ControlPoints[j][1];
        }
        Spline2d spline = Eigen::SplineFitting<Spline2d>::Interpolate(pts, 3);
        Splines.push_back(spline);
    }
}
void BoundaryCurve::setPieceWiseBelongingIfo(const Grid &grid, double physical_tol, double para_tol, int max_iter)
{
    for(const auto& spline : Splines)
    {
        CurveBelonging curveBelonging;
        curveBelonging.AdaptiveCheck(grid, spline);
        curveBelonging.PieceWiseBelonging(grid, spline, physical_tol, para_tol, max_iter);
        int n_cutcell = curveBelonging.getMultiIndices().size();
        PieceWiseBelongingIfo.push_back(curveBelonging);
    }
}