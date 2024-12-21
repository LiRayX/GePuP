#pragma once


#include "Vec.h"
#include "Grid.h"
#include "MultiIndexSet.h"

enum class type
{
    Dirichlet,
    Neumann,
};


/// @brief Coefficients for Laplacian operator
/// @param i Manhattan distance for the current cell
/// @param h grid.h
/// @return coefficients
double LaplacainCoefficients(int i, double h);
/// @brief Coefficients for Gradient operator
/// @param i Manhattan distance for the current cell
/// @param h grid.h
/// @return coefficients
double GradientCoefficients(int i, double h);

/**Note: case counts form the boundary cell, neither the cell we are computering, nor the ghost cell it self***/
double GhostCellCoefficients(int layer, int i, double h, type t);

/// @brief Ghost cell coefficients for Dirichlet boundary condition
/// @param layer Two layers of ghost cell, 1 or 2
/// @param i i=0 means the cell on the boudary, i=1 means near the boundary, i=-1 means the coffecient for boundary face
/// @return coefficients
double GhostCellCoefficients_Dirichlet(int layer, int i);
/// @brief Ghost cell coefficients for Neumann boundary condition
/// @param layer Two layers of ghost cell, 1 or 2
/// @param i i=0 means the cell on the boudary, i=1 means near the boundary, i=-1 means the coffecient for boundary face
/// @param h grid.h
/// @return coefficients
double GhostCellCoefficients_Neumann(int layer, int i, double h);
/// @brief Coputering the normal derivative on boudary face from the cell average and the boundary face averge
/// @param i 0,1,2,3 means relative Manhattan distance for the boundary face, -1 means the boundary face
/// @param h grid.h
/// @return coefficients 
double FirstNormalDerivativeCoefficients(int i, double h);
/// @brief Coputering the 2nd normal derivative on boudary face from the cell average and the boundary face averge
/// @param i 0,1,2,3 means relative Manhattan distance for the boundary face, -1 means the boundary face
/// @param h grid.h
/// @return coefficients 
double SecondNormalDerivativeCoefficients(int i, double h);

/// @brief Convert the cell average to the face average
/// @param i Manhattan distance for the computering face
/// @return coefficients
double CellAvr2FaceAvrCoefficients(int i);
double CellAvr2NormalFaceAvrCoefficients(int i, double h);


/**f(x) = \sum_j f(x_j)*L_j(x)  \sum_j f(x_j) \times\frac{\prod_{i \neq j} (x-x_i) }{\prod_{i \neq j}(x_j-x_i)} *****/
/// @brief Natural extra plotation coefficients, generated from Lagarian Polynomial Interpolation, preparation for natural outflow or extra plotation in computering Duu
/// @param i Manhattan distance for the computering cell/face
/// @return coefficients
double NaturalExtraPlotationCoefficients(int i);

double GhostCellCoefficients(int layer, int i, double h, type t)
{
    if (t == type::Dirichlet)
    {
        return GhostCellCoefficients_Dirichlet(layer, i);
    }
    else if (t == type::Neumann)
    {
        return GhostCellCoefficients_Neumann(layer, i, h);
    }
    else
    {
        std::cerr << "Not a valid type" << std::endl;
        return 0.0;
    }
}

double GhostCellCoefficients_Dirichlet(int layer, int i)
{
    if (layer == 1)
    {
        switch (i)
            {
            case 0:
                return -13.0/3;
                break;
            case 1:
                return 5.0/3;
                break;
            case 2:
                return -1.0/3;
                break;
            case -1:
                return 4.0;
                break;
            default:
                return 0.0;
                break;
    }
    }
    else if (layer == 2)
    {
       switch (i)
            {
            case 0:
                return -70.0/3;
                break;
            case 1:
                return 32.0/3;
                break;
            case 2:
                return -7.0/3;
                break;
            case -1:
                return 16.0;
                break;
            default:
                return 0.0;
                break;
            }
    }
    else
    {
        std::cerr << "Not a ghost cell" << std::endl;
        return 0.0;
    }
}

double GhostCellCoefficients_Neumann(int layer, int i, double h)
{
    if (layer == 1)
    {
        switch (i)
            {
            case 0:
                return 5.0/10;
                break;
            case 1:
                return 9.0/10;
                break;
            case 2:
                return -5.0/10;
                break;
            case 3:
                return 1.0/10;
                break;
            case -1:
                return 6.0/5 * h;
                break;
            default:
                return 0.0;
                break;
    }
    }
    else if (layer == 2)
    {
       switch (i)
            {
            case 0:
                return -75.0/10;
                break;
            case 1:
                return 145.0/10;
                break;
            case 2:
                return -75.0/10;
                break;
            case 3:
                return 15.0/10;
                break;
            case -1:
                return 6.0*h;
                break;
            default:
                return 0.0;
                break;
            }
    }
    else
    {
        std::cerr << "Not a ghost cell" << std::endl;
        return 0.0;
    }
}


double LaplacainCoefficients(int i, double h)
{
    //-1 16 -30 16 -1
    switch (i)
    {
    case 0:
        return -30.0 / (12*h * h);
        break;
    case 1:
        return 16.0 / (12*h * h);
        break;
    case 2:
        return -1.0 / (12*h * h);
        break;
    case -1:
        return 16.0 / (12*h * h);
        break;
    case -2:
        return -1.0 / (12*h * h);
        break;
    default:
        return 0.0;
        break;
    }
}

double GradientCoefficients(int i, double h)
{
    //-1 8 0 -8 1
    switch (i)
    {
    case 0:
        return 0.0;
        break;
    case 1:
        return 8.0 / (12*h);
        break;
    case 2:
        return -1.0 / (12*h);
        break;
    case -1:
        return -8.0 /(12*h);
        break;
    case -2:
        return 1.0 / (12*h);
        break;
    default:
        return 0.0;
        break;
    }
}


double FirstNormalDerivativeCoefficients(int i, double h)
{
    switch (i)
    {
        case 0:
            return -415.0/(72*h);
            break;
        case 1:
            return 161.0/(72*h);
            break;
        case 2:
            return -55.0/(72*h);
            break;
        case 3:
            return 9.0/(72*h);
            break;
        case -1:
            return 300.0/(72*h);
            break;
        default:
            return 0.0;
            break;
    }
}


double SecondNormalDerivativeCoefficients(int i, double h)
{
    switch (i)
    {
        case 0:
            return -755.0/(48*h*h);
            break;
        case 1:
            return 493.0/(48*h*h);
            break;
        case 2:
            return -191.0/(48*h*h);
            break;
        case 3:
            return 33.0/(48*h*h);
            break;
        case -1:
            return 420.0/(48*h*h);
            break;
        default:
            return 0.0;
            break;
    }

}

double CellAvr2FaceAvrCoefficients(int i)
{
    switch (i)
    {
        case -1:
            return 7.0/12;
            break;
        case 1:
            return 7.0/12;
            break;
        case -2:
            return -1.0/12;
            break;
        case 2:
            return -1.0/12;
            break;
        default:
            return 0.0;
            break;
    }
}

double CellAvr2NormalFaceAvrCoefficients(int i, double h)
{
    switch (i)
    {
        case -1:
            return -15.0/(12*h);
            break;
        case 1:
            return 15.0/(12*h);
            break;
        case -2:
            return 1.0/(12*h);
            break;
        case 2:
            return -1.0/(12*h);
            break;
        default:
            return 0.0;
            break;
    }
}


double NaturalExtraPlotationCoefficients(int i)
{
    switch (i)
    {
        case 1:
            return 5.0;
            break;
        case 2:
            return -10.0;
            break;
        case 3:
            return 10.0;
            break;
        case 4:
            return -5.0;
            break;
        case 5:
            return 1.0;
            break;
        default:
            return 0.0;
            break;
    }
}