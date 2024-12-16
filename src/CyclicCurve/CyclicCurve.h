#pragma once


#include "../../lib/Eigen/Core"
#include "../Vec.h"
#include "../Grid.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <set>

const double Pi = 4 * std::atan(1);

using ParaSet = std::set<double>;


class CyclicCurve
{
public:
    CyclicCurve(Vec _center, double radius) : center(_center), radius(radius) {}

    void setCenter(Vec _center) { center = _center; }
    void setRadius(double _radius) { radius = _radius; }

    const Vec &getCenter() const { return center; }
    double getRadius() const { return radius; }


    Vec getNormal(const Vec &point) const;
    Vec getTangent(const Vec &point) const;
    
    Vec getPoint(double theta) const;
    Vec getDer(double theta) const;
    Vec getNormal(double theta) const;
    Vec getTangent(double theta) const;

    void setIntersections(const Grid &grid);
    const ParaSet &getIntersections() const { return intersections; }

    ~CyclicCurve() = default;
protected:
    Vec center;
    double radius;
    ParaSet intersections;
};


Vec CyclicCurve::getNormal(const Vec &point) const
{
    Vec normal = normalize(point - center);
    return normal;
}

Vec CyclicCurve::getTangent(const Vec &point) const
{
    Vec tangent = counterclockwise(getNormal(point));
    return tangent;
}

Vec CyclicCurve::getPoint(double theta) const
{
    return center + Vec{std::cos(theta), std::sin(theta)} * radius;
}

Vec CyclicCurve::getDer(double theta) const
{
    return Vec{-std::sin(theta), std::cos(theta)} * radius;
}

Vec CyclicCurve::getNormal(double theta) const
{
    return getNormal(getPoint(theta));
}

Vec CyclicCurve::getTangent(double theta) const
{
    return getTangent(getPoint(theta));
}

void CyclicCurve::setIntersections(const Grid &grid)
{
    Vec left = center - Vec{radius, 0.0};
    Vec right = center + Vec{radius, 0.0};
    Vec up = center + Vec{0.0, radius};
    Vec down = center - Vec{0.0, radius};

    double h = grid.get_h();
    Vec lo = grid.lo();

    double x_origin = lo[0];
    double y_origin = lo[1];

    MultiIndex centerIndex = grid.LocateCell(center); 

    MultiIndex leftIndex = grid.LocateCell(left);
    MultiIndex rightIndex = grid.LocateCell(right);
    MultiIndex upIndex = grid.LocateCell(up);
    MultiIndex downIndex = grid.LocateCell(down);

    int i_left = leftIndex[0] + 1;
    int i_right = rightIndex[0];


    int j_down = downIndex[1] + 1;
    int j_up = upIndex[1];

    for(size_t i = i_left; i <= i_right; i++)
    {  
        double x_line = x_origin + i * h;
        double theta = std::acos((x_line - center[0]) / radius);
        intersections.insert(theta);
        intersections.insert(2*Pi-theta);
        
    }

    for(size_t j = j_down; j <= j_up; j++)
    {
        double y_line = y_origin + j * h;
        double theta = std::asin((y_line - center[1]) / radius);
        if (theta > 0)
        {
            intersections.insert(theta);
            intersections.insert(Pi-theta);
        }
        else
        {
            intersections.insert(Pi+std::fabs(theta));
            intersections.insert(2*Pi-std::fabs(theta));
        }
    }
}