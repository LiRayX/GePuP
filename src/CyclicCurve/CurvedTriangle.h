#pragma once
#include "CyclicCurve.h"
#include "../numlib.h"
#include "../Vec.h"
#include <iostream>

using ParaInterval = std::pair<double, double>;
using VecList = std::vector<Vec>;

class Triangle;
class CurvedTriangle;
class CurvedQuadrilateral;

class Triangle
{
public:
    Triangle() = default;

    Triangle(const VecList &_vertices) : vertices(_vertices) {}
    Triangle(const Vec &_v0, const Vec &_v1, const Vec &_v2) : vertices(VecList{_v0, _v1, _v2}) {}

    // void setVertices(const VecList &_vertices) { vertices = _vertices; }
    const VecList &getVertices() const { return vertices; }
    
    Vec getCentroid() const 
    {
        return (vertices[0] + vertices[1] + vertices[2]) / 3;
    }

    double getVolume() const 
    {
        return 0.5 * std::fabs(det(vertices[1] - vertices[0], vertices[2] - vertices[0]));
    }

    ~Triangle() = default;

protected:
    VecList vertices;
};


class CurvedTriangle
{
public:
    CurvedTriangle() = default;
    //constructor
    CurvedTriangle(const CyclicCurve &_cycle, const Vec &_vertex) : cycle(_cycle), vertex(_vertex) {}
    CurvedTriangle(const CyclicCurve &_cycle) : cycle(_cycle) {}
    CurvedTriangle(const CyclicCurve &_cycle, const Vec &_vertex, double _lambda_0, double _lambda_1) : cycle(_cycle), vertex(_vertex), lambda{_lambda_0, _lambda_1} {}
    CurvedTriangle(const CyclicCurve &_cycle, const Vec &_vertex, const ParaInterval &_lambda) : cycle(_cycle), vertex(_vertex), lambda(_lambda) {}

    void setVertex(Vec _vertex) { vertex = _vertex; }
    const Vec &getVertex() const { return vertex; }

    void setLambda(double _lambda_0, double _lambda_1) 
    {
        lambda.first = _lambda_0;
        lambda.second = _lambda_1;
    }

    // Lambda function
    /// @brief Jacobian of the transformation
    auto Jacobian() const;
    double Jacobian(const Vec &para) const;

    /// @brief Computering centroid of the curved triangle
    /// @return Vec of the centroid
    auto CentroidIntegral() const;
    Vec CentroidIntegral(const Vec &para) const;

    /// @brief Get the volume of the curved triangle
    /// @return double
    double getVolume() const;

    /// @brief Get the centroid of the curved triangle
    /// @return Vec of two components
    Vec getCentroid() const;
 


    /// @brief Scalar version
    /// @return 
    // auto Centroid_X() const;
    // double Centroid_X(Vec para) const;
    // auto Centroid_Y() const;
    // double Centroid_Y(Vec para) const;
    // ~CurvedTriangle() = default;
protected:
    const CyclicCurve &cycle;
    Vec vertex;
    ParaInterval lambda;
    // double lambda[2];
    // double lambda_1;
    // double lambda_2;
};


class CurvedQuadrilateral
{
public:
    CurvedQuadrilateral() = default;

    CurvedQuadrilateral(const CyclicCurve &_cycle, const Vec &_vertex_1, const Vec &_vertex_2, ParaInterval _lambda) : cycle(_cycle), vertex_1(_vertex_1), vertex_2(_vertex_2), lambda(_lambda) 
    {
        //Default counter clockwise v1 i1 i2 v2
        //curvedTriangle: v1 i1 i2
        //triangle: v1 v2 i2
        //diagnoal: v1i2, v2i1
        AdjustVertices();
        CurvedTriangle curvedTriangle(cycle, vertex_1, lambda.first, lambda.second);
        Triangle triangle(vertex_1, vertex_2, cycle.getPoint(lambda.second));
        volume = curvedTriangle.getVolume() + triangle.getVolume();
        centroid = (curvedTriangle.getCentroid() * curvedTriangle.getVolume() + triangle.getCentroid() * triangle.getVolume()) / volume;
    }

    void setVertex_1(const Vec &_vertex_1) { vertex_1 = _vertex_1; }
    void setVertex_2(const Vec &_vertex_2) { vertex_2 = _vertex_2; }
    void setLambda(ParaInterval _lambda) { lambda = _lambda; }

    void setLambda(double _lambda_0, double _lambda_1) 
    {
        lambda.first = _lambda_0;
        lambda.second = _lambda_1;
    }

    const Vec &getVertex_1() const { return vertex_1; }
    const Vec &getVertex_2() const { return vertex_2; }
    const ParaInterval &getLambda() const { return lambda; }
    /// @brief Adjust the order of the vertices to form a quadrilateral
    void AdjustVertices();
    /// @brief Get the volume of the curved triangle
    /// @return double
    double getVolume() const { return volume; }

    /// @brief Get the centroid of the curved triangle
    /// @return Vec of two components
    Vec getCentroid() const { return centroid; }

protected:
    const CyclicCurve &cycle;
    Vec vertex_1;
    Vec vertex_2;
    ParaInterval lambda;
    double volume;
    Vec centroid;
};




void CurvedQuadrilateral::AdjustVertices()
{
    Vec intersection_1 = cycle.getPoint(lambda.first);
    Vec intersection_2 = cycle.getPoint(lambda.second);
    //Range of argument is [-pi,pi]
    double argu_v1v2 = argument(vertex_2 - vertex_1);
    double argu_v1i1 = argument(intersection_1 - vertex_1);
    double argu_v1i2 = argument(intersection_2 - vertex_1);

    std::vector<double> arguments = {argu_v1v2, argu_v1i1, argu_v1i2};
    std::sort(arguments.begin(), arguments.end());
    //According to relative position to the line of -x-axis, where argument = -pi or pi
    if(arguments[2] - arguments[0] > Pi)
    {
        if(arguments[0]<0 && arguments[1]<0 && arguments[2]>0)
        {
            if(argu_v1i1 == arguments[0])
            {
                swap(vertex_1, vertex_2);
            }
        }
        else if(arguments[0]<0 && arguments[1]>0 && arguments[2]>0)
        {
            if(argu_v1i1 == arguments[2])
            {
                swap(vertex_1, vertex_2);
            }
        }
    }
    //If v1i1 is the diagonal
    else
    {
        if(argu_v1i1 == arguments[1])
            std::cout << "v1i1 is the diagonal" << std::endl;
            swap(vertex_1, vertex_2);
    }
}






auto CurvedTriangle::Jacobian() const 
{
    return [this](const Vec& para) -> double 
    {
        return this->Jacobian(para);
    };
}
double CurvedTriangle::Jacobian(const Vec &para) const
{
    // X = (1-\theta)C(\lambda) + \theta V
    //C(\lambda)
    Vec point = cycle.getPoint(para[0]);
    //C'(\lambda)
    Vec der = cycle.getDer(para[0]) * (1 - para[1]);
    double jacobian = std::fabs(det(der, vertex - point));
    return jacobian;
}

auto CurvedTriangle::CentroidIntegral() const
{
    return [this](const Vec& para) -> Vec 
    {
        return this->CentroidIntegral(para);
    };
}
Vec CurvedTriangle::CentroidIntegral(const Vec &para) const
{
    //X = (1-\theta)C(\lambda) + \theta * V
    Vec coord = cycle.getPoint(para[0])*(1 - para[1]) + vertex*para[1];  
    return coord * Jacobian(para);
}

double CurvedTriangle::getVolume() const
{
    return quad2D_scalar(Jacobian(), lambda.first, lambda.second, 0, 1);
}

Vec CurvedTriangle::getCentroid() const
{
    Vec centroid = quad2D_vector(CentroidIntegral(), lambda.first, lambda.second, 0, 1);
    return centroid / getVolume();
}


// auto CurvedTriangle::Centroid_X() const
// {
//     return [this](const Vec& para) -> double 
//     {
//         return this->Centroid_X(para);
//     };
// }
// double CurvedTriangle::Centroid_X(Vec para) const
// {
//     Vec point = cycle.getPoint(para[0]);
//     Vec der = cycle.getDer(para[0]) * (1 - para[1]);
//     double jacobian = std::fabs(det(der, vertex - point));
//     return jacobian * point[0];
// }

// auto CurvedTriangle::Centroid_Y() const
// {
//     return [this](const Vec& para) -> double 
//     {
//         return this->Centroid_Y(para);
//     };
// }
// double CurvedTriangle::Centroid_Y(Vec para) const
// {
//     Vec point = cycle.getPoint(para[0]);
//     Vec der = cycle.getDer(para[0]) * (1 - para[1]);
//     double jacobian = std::fabs(det(der, vertex - point));
//     return jacobian * point[1];
// }