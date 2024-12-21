#pragma once

#include "Vec.h"
#include "Grid.h"
#include "MultiIndexSet.h"
#include "Data.h"



ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, int max_iter);

ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, double tol);

ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, int max_iter)
{
    Eigen::VectorXd D = A.diagonal();
    Eigen::VectorXd D_inv = D.cwiseInverse();
    ScalarData x = x0;
    for(int k=0; k<max_iter; k++)
    {
        Eigen::VectorXd residual = (b - Multiply(A,x)).get_data();
        x.get_data() += weight * D_inv.cwiseProduct(residual);
    }
    return x;
}

ScalarData WeightedJacobi(const DiffMatrix &A, const ScalarData &b, const ScalarData &x0, double weight, double tol)
{
    Eigen::VectorXd D = A.diagonal();
    Eigen::VectorXd D_inv = D.cwiseInverse();
    ScalarData x = x0;
    Eigen::VectorXd residual = (b - Multiply(A,x)).get_data();
    while(residual.norm() > tol)
    {
        x.get_data() += weight * D_inv.cwiseProduct(residual);
        residual = (b - Multiply(A,x)).get_data();
    }
    return x;
}
