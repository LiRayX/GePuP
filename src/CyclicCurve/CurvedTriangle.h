#pragma once
#include "CyclicCurve.h"
#include "../numlib.h"

using ParaInterval = std::pair<double, double>;

class CurvedTriangle
{
public:
    //constructor
    CurvedTriangle(const CyclicCurve &_cycle, Vec _vertex) : cycle(_cycle), vertex(_vertex) {}
    CurvedTriangle(const CyclicCurve &_cycle) : cycle(_cycle) {}
    CurvedTriangle(const CyclicCurve &_cycle, Vec _vertex, double _lambda_0, double _lambda_1) : cycle(_cycle), vertex(_vertex), lambda{_lambda_0, _lambda_1} {}
    CurvedTriangle(const CyclicCurve &_cycle, Vec _vertex, ParaInterval _lambda) : cycle(_cycle), vertex(_vertex), lambda(_lambda) {}

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
    double Jacobian(Vec para) const;

    /// @brief Computering centroid of the curved triangle
    /// @return Vec of the centroid
    auto CentroidIntegral() const;
    Vec CentroidIntegral(Vec para) const;

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
class Triangle
{
public:
    Triangle() = default;

    Triangle(const VecList &_vertices) : vertices(_vertices) {}
    Triangle(const Vec &_v0, const Vec &_v1, const Vec &_v2) : vertices(VecList{_v0, _v1, _v2}) {}

    // void setVertices(const VecList &_vertices) { vertices = _vertices; }
    const VecList &getVertices() const { return vertices; }
    
    Vec Centroid() const 
    {
        return (vertices[0] + vertices[1] + vertices[2]) / 3;
    }

    double Volume() const 
    {
        return 0.5 * std::fabs(det(vertices[1] - vertices[0], vertices[2] - vertices[0]));
    }

    ~Triangle() = default;

protected:
    const VecList &vertices;
};



auto CurvedTriangle::Jacobian() const 
{
    return [this](const Vec& para) -> double 
    {
        return this->Jacobian(para);
    };
}
double CurvedTriangle::Jacobian(Vec para) const
{
    Vec point = cycle.getPoint(para[0]);
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
Vec CurvedTriangle::CentroidIntegral(Vec para) const
{
    Vec point = cycle.getPoint(para[0]);
    Vec der = cycle.getDer(para[0]) * (1 - para[1]);
    double jacobian = std::fabs(det(der, vertex - point));
    return point * jacobian;
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