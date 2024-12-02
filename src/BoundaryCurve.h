#pragma once


#include "../lib/Eigen/Core"
#include "../lib/unsupported/Eigen/Splines"
#include "Vec.h"
#include "Grid.h"
#include <iostream>
#include <vector>

using VecList = std::vector<Vec>;
using Spline2d = Eigen::Spline<double, 2, 3>;
using SplineList = std::vector<Spline2d>;

class BoundaryCurve
{
public:
    BoundaryCurve() = default;

    BoundaryCurve(const VecList& controlPoints, const std::vector<int>& endPoints)
        : ControlPoints(controlPoints), IndexEndPoints(endPoints) {}

    VecList GetControlPoints() const { return ControlPoints; }
    std::vector<int> GetEndPoints() const { return IndexEndPoints; }

    SplineList GetSplines() const { return Splines; }

    
    Spline2d GlobalInterpolate() const;
    void LocalInterpolate();




protected:
    VecList ControlPoints;
    std::vector<int> IndexEndPoints;
    SplineList Splines;
};



Spline2d BoundaryCurve::GlobalInterpolate() const
{
    Eigen::MatrixXd pts(2, ControlPoints.size());
    for (int i = 0; i < ControlPoints.size(); ++i)
    {
        pts(0, i) = ControlPoints[i][0];
        pts(1, i) = ControlPoints[i][1];
    }

    return Eigen::SplineFitting<Spline2d>::Interpolate(pts, 3);
}

void BoundaryCurve::LocalInterpolate()
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