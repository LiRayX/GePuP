// #define EIGEN_MALLOC_ALREADY_ALIGNED 0
#include "../lib/Eigen/Core"
#include "../lib/unsupported/Eigen/Splines"
#include "../src/Vec.h"
#include "../src/BoundaryCurve.h"
#include <iostream>


using VecList = std::vector<Vec>;
using Spline2d = Eigen::Spline<double, 2, 3>;
using SplineList = std::vector<Spline2d>;

int main()
{
    VecList points;
    points.push_back(Vec{0, 0});
    points.push_back(Vec{1, 2});
    points.push_back(Vec{3, 1});
    points.push_back(Vec{4, 0});
    points.push_back(Vec{5, 0});
    points.push_back(Vec{6, 2});
    points.push_back(Vec{6, 0});

    std::vector<int> endPoints = {0, 3, 6};

    BoundaryCurve boundaryCurve(points, endPoints);
    boundaryCurve.LocalInterpolate();

    SplineList splines = boundaryCurve.GetSplines();
    for (int i = 0; i < splines.size(); ++i)
    {
        for (double u = 0; u <= 1; u += 0.1)
        {
            Eigen::Vector2d point = splines[i](u);
            std::cout << "u = " << u << ", point = (" << point(0) << ", " << point(1) << ")" << std::endl;
        }
    }

    return 0;
}