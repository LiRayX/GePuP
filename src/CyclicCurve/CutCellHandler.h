#pragma once

#include "../Vec.h"
#include "../Grid.h"
#include "../Function.h"
#include "../MultiIndexSet.h"
#include "../numlib.h"

#include "CellDivision.h"
#include "CyclicCurve.h"
#include <vector>

using VecList = std::vector<Vec>;
using CutCellMapping = std::unordered_map<MultiIndex, std::set<double>>;

void swap(double &a, double &b);

class CutCellHandler
{
public:
    int NumOutsideCorner(const MultiIndex &index, const Grid &grid, const CyclicCurve &boundaryCurve);
    VecList getOutsideCorners(const MultiIndex &index, const Grid &grid, const CyclicCurve &boundaryCurve);
    double Volume(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);


};


int CutCellHandler::NumOutsideCorner(const MultiIndex &index, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    int num = 0;
    VecList corners = {grid(index[0], index[1]), grid(index[0] + 1, index[1]), grid(index[0] + 1, index[1] + 1), grid(index[0], index[1] + 1)};
    for (const auto &corner : corners)
    {
        if (norm(corner-boundaryCurve.getCenter()) > boundaryCurve.getRadius())
        {
            num++;
        }
    }
    return num;
}
VecList CutCellHandler::getOutsideCorners(const MultiIndex &index, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    VecList outside_corners;
    VecList corners = {grid(index[0], index[1]), grid(index[0] + 1, index[1]), grid(index[0] + 1, index[1] + 1), grid(index[0], index[1] + 1)};
    for (const auto &corner : corners)
    {
        if (norm(corner-boundaryCurve.getCenter()) > boundaryCurve.getRadius())
        {
            outside_corners.push_back(corner);
        }
    }
    return outside_corners;
}



double CutCellHandler::Volume(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    //Init the volume
    double volume = 0;
    //Basic information of grid
    double h = grid.get_h();
    //Get the outside corners
    VecList outside_corners = getOutsideCorners(index, grid, boundaryCurve);
    int num_outside_corner = outside_corners.size();
    //For the last cut-cell, it's para range span the x-axis
    double theta_0 = *para.begin();
    double theta_1 = *para.rbegin();
    if (theta_1 - theta_0 > Pi)
    {
        theta_1 -= 2 * Pi;
    }
    //Ensure theta_0 < theta_1
    swap(theta_0, theta_1);
    //According to the number of outside corners, we can determine the volume of the cut-cell
    switch (num_outside_corner)
    {
        //All inside the circle, this won't happen since the cycle is convex
    case 0:
        break;
        //One outside the circle, a curved triangle
    case 1:
        break;
        //Two outside the circle, a curved quadrilateral
    case 2:
        volume = 0.5 * h * h;
        break;
        //Complement of case 1
    case 3:
        volume = 0.25 * h * h;
        break;
        //Cut the same edge twice
    case 4:
        volume = 0;
        break;
    default:
        break;
    }
    return volume;
}

void swap(double& a, double& b)
{
    double temp = a;
    a = b;
    b = temp;
}

Vec ParaCurve(const Vec &center, double radius, double lambda)
{
    return center + Vec{std::cos(lambda), std::sin(lambda)}*radius;
}
