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
/// @brief Adjust the parameter interval
ParaInterval AdjustParaInterval(const ParaSet &para);
/// @brief Get the direction to find the vertex of the triangle
Vec getDirection(const Vec &intersection, const Vec &center);

class CutCellHandler
{
public:
    CutCellHandler() = default;
    void ClassifyCorners(const MultiIndex &index, const Grid &grid, const CyclicCurve &boundaryCurve);
    void Handler(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);

    void OneInsideCorner(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);

    void TwoInsideCorner2(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);
    void TwoInsideCorner4(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);
    void ThreeInsideCorner(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);

    void FourInsideCorner(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve);


    void setVolume(double _volume) { volume = _volume; }
    double getVolume() const { return volume; }
    void setCentroid(Vec _centroid) { centroid = _centroid; }
    Vec getCentroid() const { return centroid; }


protected:
    double volume;
    Vec centroid;
    VecList inside_corners;
    VecList outside_corners;
};



void CutCellHandler::ClassifyCorners(const MultiIndex &index, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    VecList corners = {grid(index[0], index[1]), grid(index[0] + 1, index[1]), grid(index[0] + 1, index[1] + 1), grid(index[0], index[1] + 1)};
    for (const auto &corner : corners)
    {
        if ((norm(corner-boundaryCurve.getCenter()) - boundaryCurve.getRadius()) < 0)
        {
            inside_corners.push_back(corner);
        }
        else
        {
            outside_corners.push_back(corner);
        }
    }
}



void CutCellHandler::Handler(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    //If only one intersection point, then it is a whole cell
    if (para.size() == 1)
    {
        volume = grid.get_cell_volume();
        centroid = grid.center(index);
    }
    //Init the volume
    double volume = 0;
    //Basic information of grid
    double h = grid.get_h();
    //Get the outside corners
    //The number of outside corners
    int num_inside_corner = inside_corners.size();
    //According to the number of outside corners, we can determine the volume of the cut-cell
    switch (num_inside_corner)
    {
        //All outside the circle, if this happens, then a neithbor is cut by the circle fourth times
    case 0:
        break;
        //One inside the circle, a curved triangle
    case 1:
        OneInsideCorner(index, para, grid, boundaryCurve);
        break;
        //Two outside the circle, a curved quadrilateral
    case 2:
        if (para.size() == 2)
            TwoInsideCorner2(index, para, grid, boundaryCurve);
        else if (para.size() == 4)
            TwoInsideCorner4(index, para, grid, boundaryCurve);
        break;
        //Three inside the circle, a curved pentagon
    case 3:
        ThreeInsideCorner(index, para, grid, boundaryCurve);
        break;
        //All inside the circle, , this won't happen since the cycle is convex
    case 4:
        FourInsideCorner(index, para, grid, boundaryCurve);
        break;
    default:
        break;
    }
}


void CutCellHandler::OneInsideCorner(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    ParaInterval lambda = AdjustParaInterval(para);
    /**********************************Computering the Curved Quadrilateral*********************************/
    //In this case, the convex curved quadrilateral is divided into a curved triangle and a triangle
    //And the alive region is the Complement of the whole cell
    Vec corner = *inside_corners.begin();
    //Construct the curved triangle
    CurvedTriangle curvedTriangle(boundaryCurve, corner, lambda);
    //Get the volume of the alive region
    volume = grid.get_cell_volume() - curvedTriangle.getVolume();
    //Get the centroid integral of the alive region
    centroid = grid.center(index)*grid.get_cell_volume() - curvedTriangle.getCentroid()*curvedTriangle.getVolume();
    //Get the centroid of the alive region
    centroid = centroid/volume;
}



void CutCellHandler::TwoInsideCorner2(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{

    ParaInterval lambda = AdjustParaInterval(para);
    /**********************************Computering the Curved Quadrilateral*********************************/
    //In this case, the convex curved quadrilateral is divided into a curved triangle and a triangle
    //And the alive region is the Complement of the whole cell
    Vec vertex = *inside_corners.begin();
    Vec vertex_t = *inside_corners.rbegin();
    //Set the vertex of the curved triangle, each inside corner is OK
    //Construct the curved triangle
    CurvedTriangle curvedTriangle(boundaryCurve, vertex, lambda);
    //Set the vertices of the triangle
        //Get the intersection points
    Vec intersection = boundaryCurve.getPoint(lambda.first);
    //If the intersection point is on the same face as the vertex of the curved triangle
    //Then the left point will be a vertex of the triangle
    Vec diff = intersection - vertex;
    double tol = 1e-12;
    if (std::fabs(diff[0])<tol || std::fabs(diff[1])<tol)
    {
        intersection = boundaryCurve.getPoint(lambda.second);
    }
    Triangle triangle(vertex, vertex_t, intersection);

    volume = grid.get_cell_volume() - triangle.Volume() - curvedTriangle.getVolume();
    centroid = grid.center(index)*grid.get_cell_volume() - triangle.Centroid()*triangle.Volume() - curvedTriangle.getCentroid()*curvedTriangle.getVolume();
        //
    //Get the centroid of the alive region
    centroid = centroid/volume;
}

void CutCellHandler::TwoInsideCorner4(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    auto it = para.begin();
    ParaInterval lambda_interval_1 = std::make_pair(*it, *(it++));
    ParaInterval lambda_interval_2 = std::make_pair(*(it++), *(it++));

    ParaInterval lambda_interval_1 = AdjustParaInterval(lambda_interval_1);
    ParaInterval lambda_interval_2 = AdjustParaInterval(lambda_interval_2);
    /**********************************Computering the Curved Quadrilateral*********************************/
    //In this case, the convex curved quadrilateral is divided into a curved triangle and a triangle
    //And the alive region is the Complement of the whole cell
    Vec vertex = *inside_corners.begin();
    Vec vertex_t = *inside_corners.rbegin();
    //Set the vertex of the curved triangle, each inside corner is OK
    //Construct the curved triangle
    CurvedTriangle curvedTriangle_1(boundaryCurve, vertex, lambda_interval_1);
    CurvedTriangle curvedTriangle_2(boundaryCurve, vertex, lambda_interval_2);
    //Set the vertices of the triangle
        //Get the intersection points
    Vec intersection_mid_1 = boundaryCurve.getPoint(lambda_interval_1.second);
    Vec intersection_mid_2 = boundaryCurve.getPoint(lambda_interval_2.first);
    //Find the corner away from the outside corner, this corner will be the vertex of triangle
    Vec intersection = boundaryCurve.getPoint(lambda_interval_1.first);
    //If the intersection point is on the same face as the vertex of the curved triangle
    //Then the left point will be a vertex of the triangle
    Vec diff = intersection - vertex;
    double tol = 1e-12;
    if (std::fabs(diff[0])<tol || std::fabs(diff[1])<tol)
    {
        intersection = boundaryCurve.getPoint(lambda_interval_2.second);
    }
    //Construct the triangle
    Triangle triangle_1(vertex, intersection_mid_1, intersection_mid_2);
    Triangle triangle_2(vertex, vertex_t, intersection);

    volume = grid.get_cell_volume() - triangle_1.Volume() - triangle_2.Volume() - curvedTriangle_1.getVolume() - curvedTriangle_2.getVolume();
    centroid = grid.center(index)*grid.get_cell_volume() - triangle_1.Centroid()*triangle_1.Volume() - triangle_2.Centroid()*triangle_2.Volume() - curvedTriangle_1.getCentroid()*curvedTriangle_1.getVolume() - curvedTriangle_2.getCentroid()*curvedTriangle_2.getVolume();
    //Get the centroid of the alive region
    centroid = centroid/volume;
}


void CutCellHandler::ThreeInsideCorner(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    ParaInterval lambda = AdjustParaInterval(para);

    /**********************************Computering the Curved Pentagon*********************************/
    //In this case, the curved pentagon is divided into a curved triangle and two triangles
    //And the alive region is the Complement of the whole cell
        //step-1: Get the vertex from the outside corner along the diagonal direction AND Construct the curved triangle
    Vec outside_corner = *outside_corners.begin();
    Vec vertex = grid.center(index)*2 - outside_corner;
    CurvedTriangle curvedTriangle(boundaryCurve, vertex, lambda);

        //step-2: Specify the two triangles
    //Get Intersection points
    Vec intersection_1 = boundaryCurve.getPoint(lambda.first);
    Vec intersection_2 = boundaryCurve.getPoint(lambda.second);

    //Find the corner away from the outside corner, this corner will be the vertex of triangle
    Vec vertex_1 = outside_corner + getDirection(intersection_1, outside_corner) * grid.get_h();
    Vec vertex_2 = outside_corner + getDirection(intersection_2, outside_corner) * grid.get_h();
    //Construct the triangle
    Triangle triangle_1(vertex_1, vertex, intersection_1);
    Triangle triangle_2(vertex_2, vertex, intersection_2);


        //step-3: Get the volume of the alive region
    //The volume of the alive region
    volume = grid.get_cell_volume() - curvedTriangle.getVolume() - triangle_1.Volume() - triangle_2.Volume();
    //Get the centroid integral of the alive region
    centroid = grid.center(index)*grid.get_cell_volume() - curvedTriangle.getCentroid()*curvedTriangle.getVolume() - triangle_1.Centroid()*triangle_1.Volume() - triangle_2.Centroid()*triangle_2.Volume();
    //Get the centroid of the alive region
    centroid = centroid/volume;
}


void CutCellHandler::FourInsideCorner(const MultiIndex &index, const ParaSet &para, const Grid &grid, const CyclicCurve &boundaryCurve)
{
    ParaInterval lambda = AdjustParaInterval(para);
    Vec center_cycle = boundaryCurve.getCenter();

    /**********************************Computering the Arc*********************************/
    //In this case, only a edge is cut by the circle

    //Get the intersection points
    Vec intersection_1 = boundaryCurve.getPoint(lambda.first);
    Vec intersection_2 = boundaryCurve.getPoint(lambda.second);
    //Construct the curved triangle with the center of the circle
    CurvedTriangle curvedTriangle(boundaryCurve, center_cycle, lambda);
    Triangle triangle  = Triangle(intersection_1, intersection_2, center_cycle);

    volume = grid.get_cell_volume() - curvedTriangle.getVolume() + triangle.Volume();
    centroid = grid.center(index)*grid.get_cell_volume() - curvedTriangle.getCentroid()*curvedTriangle.getVolume() + triangle.Centroid()*triangle.Volume();
    centroid = centroid/volume;

}   


void swap(double &a, double &b)
{
    double temp = a;
    a = b;
    b = temp;
}


ParaInterval AdjustParaInterval(const ParaInterval &para)
{
    double para_1 = para.first;
    double para_2 = para.second;
    if (para_2 - para_1 > Pi)
    {
        para_2 -= 2 * Pi;
    }
    swap(para_1, para_2);
    return std::make_pair(para_1, para_2);
}

ParaInterval AdjustParaInterval(const ParaSet &para)
{
    assert(para.size() == 2);
    ParaInterval para_interval;
    double para_1 = *para.begin();
    double para_2 = *para.rbegin();
    if (para_2 - para_1 > Pi)
    {
        para_2 -= 2 * Pi;
    }
    swap(para_1, para_2);
    para_interval.first = para_1;
    para_interval.second = para_2;
    return para_interval;
}


Vec getDirection(const Vec &intersection, const Vec &center)
{
    Vec diff = intersection - center;
    double tol = 1e-12;
    if (std::fabs(diff[0])<tol) 
    {
        double result = (diff[1] > 0) ? 1.0 : -1.0;
        return Vec{0, result};
    }
    else if (std::fabs(diff[1])<tol)
    {
        double result = (diff[0] > 0) ? 1.0 : -1.0;
        return Vec{result, 0};
    }
    else
    {
        std::cerr << "The intersection point is not on the diagonal direction" << std::endl;
        return Vec{0, 0};
    }
}




