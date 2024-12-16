#pragma once
#include "CyclicCurve.h"

class CurvedTriangle
{
public:
    // CurvedTriangle(const CyclicCurve &_cycle, Vec _vertex, double _lambda_1, double _lambda_2) : cycle(_cycle), vertex(_vertex), lambda_1(_lambda_1), lambda_2(_lambda_2) {}

    CurvedTriangle(const CyclicCurve &_cycle, Vec _vertex) : cycle(_cycle), vertex(_vertex) {}
    CurvedTriangle(const CyclicCurve &_cycle) : cycle(_cycle) {}

    void setVertex(Vec _vertex) { vertex = _vertex; }
    // void setLambda(double _lambda_1, double _lambda_2) { lambda_1 = _lambda_1; lambda_2 = _lambda_2; }
    const Vec &getVertex() const { return vertex; }
    auto Jacobian() const 
    {
        return [this](const Vec& para) -> double 
        {
            return this->Jacobian(para);
        };
    }
    double Jacobian(Vec para) const
    {
        Vec point = cycle.getPoint(para[0]);
        Vec der = cycle.getDer(para[0]) * (1 - para[1]);
        double jacobian = std::fabs(det(der, vertex - point));
        return jacobian;
    }
    double Centroid_X() const;
    ~CurvedTriangle() = default;
protected:
    const CyclicCurve &cycle;
    Vec vertex;
    // double lambda_1;
    // double lambda_2;
};