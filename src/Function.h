#pragma once

#include "../lib/Eigen/Dense"
#include <cmath>
#include "Vec.h"

using HessianMatrix = Eigen::Matrix2d;
using GradientMatrix = Eigen::Matrix2d;

class ScalarFunction
{
public:
    virtual double operator()(Vec _v) const = 0;
    // ((\partial f)/(\partial x),(\partial f)/(\partial y)) 
    virtual Vec der(Vec _v) const { return Vec{0, 0}; }
    // ((\partial^2 f)/(\partial x^2),(\partial^2 f)/(\partial y^2))
    virtual Vec der2(Vec _v) const { return Vec{0, 0}; }

    virtual ~ScalarFunction() = default; 

};

class VecFunction
{
public:
    virtual Vec operator()(Vec _v) const = 0;
    //[(\partial u_x)/(\partial x),(\partial u_x)/(\partial y)
    //(\partial u_y)/(\partial x),(\partial u_y)/(\partial y)]
    virtual GradientMatrix der(Vec _v) const { return GradientMatrix::Zero(2, 2); }
    //[(\partial^2 u_x)/(\partial x^2),(\partial^2 u_x)/(\partial y^2)
    //(\partial^2 u_y)/(\partial x^2),(\partial^2 u_y)/(\partial y^2)] 
    virtual HessianMatrix der2(Vec _v) const { return HessianMatrix::Zero(2, 2); }
    virtual ~VecFunction() = default; 
};

class PressureFeild : public ScalarFunction
{
public:
    double operator()(Vec _v) const override
    {
        return 0;
    }
    Vec der(Vec _v) const override
    {
        return Vec{0, 0};
    }
    Vec der2(Vec _v) const override
    {
        return Vec{0, 0};
    }
};
class VelocityFeild : public VecFunction
{
public:
    Vec operator()(Vec _v) const override
    {
        return Vec{0, 0};
    }
    GradientMatrix der(Vec _v) const override
    {
        GradientMatrix der(2, 2);
        return der;
    }
    HessianMatrix der2(Vec _v) const override
    {
        HessianMatrix der(2, 2);
        return der;
    }
};






// enum class Func
// {
//   pressure,
//   qforp,
//   velocity,
//   gravity
// };

// class ScalarFunction;

// double p(Vec _v);
// double q(Vec _v);

// Vec p_der(Vec _v);
// Vec q_der(Vec _v);

// Vec p_der2(Vec _v);
// Vec q_der2(Vec _v);

// class VecFunction;

// Vec u(Vec _v);
// Vec g(Vec _v);

// GradientMatrix u_der(Vec _v);
// GradientMatrix g_der(Vec _v);

// HessianMatrix u_der2(Vec _v);
// HessianMatrix g_der2(Vec _v);

// // f:R^2->R
// class ScalarFunction
// {
//   friend double p(Vec _v);
//   friend double q(Vec _v);
//   friend Vec p_der(Vec _v);
//   friend Vec q_der(Vec _v);
//   friend Vec p_der2(Vec _v);
//   friend Vec q_der2(Vec _v);

// public:
//   // Constructors
//   ScalarFunction() = default;

//   ScalarFunction(const Func _func) : f(_func) {}

//   // accessor
//   double operator()(Vec _v) const
//   {
//     if (f == Func::pressure)
//       return p(_v);
//     else if (f == Func::qforp)
//       return q(_v);
//   }
//   // ((\partial f)/(\partial x),(\partial f)/(\partial y))
//   Vec der(Vec _v) const
//   {
//     if (f == Func::pressure)
//       return p_der(_v);
//     else if (f == Func::qforp)
//       return q_der(_v);
//   }

//   // ((\partial^2 f)/(\partial x^2),(\partial^2 f)/(\partial y^2))
//   Vec der2(Vec _v) const
//   {
//     if (f == Func::pressure)
//       return p_der2(_v);
//     else
//       return q_der2(_v);
//   }

// protected:
//   Func f;
// };

// // f:R^2->R^2
// class VecFunction
// {
//   friend Vec u(Vec _v);
//   friend Vec g(Vec _v);
//   friend GradientMatrix u_der(Vec _v);
//   friend GradientMatrix g_der(Vec _v);
//   friend HessianMatrix u_der2(Vec _v);
//   friend HessianMatrix g_der2(Vec _v);

// public:
//   // Constructors
//   VecFunction() = default;

//   VecFunction(const Func _func) : f(_func) {}

//   // accessor
//   Vec operator()(Vec _v) const
//   {
//     if (f == Func::velocity)
//       return u(_v);
//     else if (f == Func::gravity)
//       return g(_v);
//   }

//   //[(\partial u_x)/(\partial x),(\partial u_x)/(\partial y)
//   //(\partial u_y)/(\partial x),(\partial u_y)/(\partial y)]
//   GradientMatrix der(Vec _v) const
//   {
//     if (f == Func::velocity)
//       return u_der(_v);
//     else if (f == Func::gravity)
//       return g_der(_v);
//   }

//   //[(\partial^2 u_x)/(\partial x^2),(\partial^2 u_x)/(\partial y^2)
//   //(\partial^2 u_y)/(\partial x^2),(\partial^2 u_y)/(\partial y^2)]
//  HessianMatrix der2(Vec _v) const
//   {
//     if (f == Func::velocity)
//       return u_der2(_v);
//     else if (f == Func::gravity)
//       return g_der2(_v);
//   }

// protected:
//   Func f;
// };
// double p(Vec _v)
// {
//   // Change the function here
//   return 1;
// };

// double q(Vec _v)
// {
//   // Change the function here
//   return 0;
// };

// Vec p_der(Vec _v)
// {
//   // Change the function here
//   return Vec{0, 0};
// };

// Vec q_der(Vec _v)
// {
//   // Change the function here
//   return Vec{0, 0};
// };

// Vec p_der2(Vec _v)
// {
//   // Change the function here
//   Vec re = _v * 2;
//   return re;
// };

// Vec q_der2(Vec _v)
// {
//   // Change the function here
//   return Vec{0, 0};
// };

// Vec u(Vec _v)
// {
//   // Change the function here
//   return 0;
// };

// Vec g(Vec _v)
// {
//   // Change the function here
//   return 0;
// };

// GradientMatrix u_der(Vec _v)
// {
//   // Change the function here
//   HessianMatrix der(2, 2);
//   return der;
// };

// GradientMatrix g_der(Vec _v)
// {
//   // Change the function here
//   HessianMatrix der(2, 2);
//   return der;
// };

// HessianMatrix u_der2(Vec _v)
// {
//   // Change the function here
//   HessianMatrix der(2, 2);
//   der<<0,0,
//       0,0;
//   return der;
// };

// HessianMatrix g_der2(Vec _v)
// {
//   // Change the function here
//   HessianMatrix der(2, 2);
//   return der;
// };