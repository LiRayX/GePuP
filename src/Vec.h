#pragma once

#include <iostream>
#include <cmath>
#include <initializer_list>

class Vec
{
public:
  // constructors
  Vec(const double &t = double())
  {
    for (int i = 0; i < 2; i++)
      coord[i] = t;
  }

  Vec(std::initializer_list<double> l)
  {
    auto j = l.begin();
    for (int i = 0; i < 2; i++)
      coord[i] = *j++;
  }

  // The conversion constructor
  Vec(const Vec &rhs)
  {
    for (int i = 0; i < 2; i++)
      coord[i] = rhs[i];
  }

  Vec &operator=(const Vec &rhs)
  {
    for (int i = 0; i < 2; i++)
      coord[i] = rhs[i];
    return *this;
  }
  // The unit vector
  static Vec unit(int D)
  {
    Vec u;
    u[D] = static_cast<double>(1);
    return u;
  }

  // accessors
  double &operator[](int d) { return coord[d]; }

  const double &operator[](int d) const { return coord[d]; }

  const double *data() const { return &coord[0]; }

public:
#define ELMWISE_BINARY_OP(OpName, Op) \
  auto OpName(const Vec &rhs) const   \
  {                                   \
    Vec res;                          \
    for (int i = 0; i < 2; i++)       \
      res[i] = coord[i] Op rhs[i];    \
    return res;                       \
  }

  ELMWISE_BINARY_OP(operator+, +)
  ELMWISE_BINARY_OP(operator-, -)
  ELMWISE_BINARY_OP(operator*, *)
  ELMWISE_BINARY_OP(operator/, /)
#undef ELMWISE_BINARY_OP

#define RIGHT_BROADCAST(OpName, Op)    \
  auto OpName(const double &rhs) const \
  {                                    \
    Vec res;                           \
    for (int i = 0; i < 2; i++)        \
      res[i] = coord[i] Op rhs;        \
    return res;                        \
  }

  RIGHT_BROADCAST(operator+, +)
  RIGHT_BROADCAST(operator-, -)
  RIGHT_BROADCAST(operator*, *)
  RIGHT_BROADCAST(operator/, /)
#undef RIGHT_BROADCAST

  Vec operator-() const
  {
    Vec res;
    for (int i = 0; i < 2; i++)
      res[i] = -coord[i];
    return res;
  }

  friend std::ostream &operator<<(std::ostream &os, const Vec &_v)
  {
    os << "(" << _v[0];
    for (int i = 1; i < 2; i++)
      os << "," << _v[i];
    os << ")";
    return os;
  }

protected:
  double coord[2];
};

// Operators of Vec

//====================================================
// binary operators

inline double dot(const Vec &lhs, const Vec &rhs)
{
  double res = 0;
  for (int i = 0; i < 2; i++)
    res += lhs[i] * rhs[i];
  return res;
}

//====================================================
// unary operators

inline Vec abs(const Vec &v)
{
  Vec res;
  for (int i = 0; i < 2; i++)
    res[i] = std::abs(v[i]);
  return res;
}

inline Vec ceiling(const Vec &v)
{
  Vec res;
  for (int i = 0; i < 2; i++)
    res[i] = std::ceil(v[i]);
  return res;
}

inline Vec floor(const Vec &v)
{
  Vec res;
  for (int i = 0; i < 2; i++)
    res[i] = std::floor(v[i]);
  return res;
}

inline Vec sign(const Vec &v)
{
  Vec res;
  for (int i = 0; i < 2; i++)
    res[i] = (v[i] > 0) ? 1 : ((v[i] < 0) ? -1 : 0);
  return res;
}

inline double norm(const Vec &v, int nt = 2)
{
  if (nt == 2)
    return sqrt(dot(v, v));
  else if (nt == 1)
  {
    double res = 0;
    for (int i = 0; i < 2; i++)
      res += std::abs(v[i]);
    return res;
  }
  else if (nt == 0)
  {
    double res = std::abs(v[0]);
    for (int i = 1; i < 2; i++)
      res = std::max(std::abs(v[i]), res);
    return res;
  }
  else
  {
    std::cerr << "Error: norm type not supported!" << std::endl;
    return 0;
  }
}

inline Vec normalize(const Vec &v)
{
  double l = norm(v, 2);
  return v / l;
}

inline double sum(const Vec &v)
{
  double res = v[0];
  for (int i = 1; i < 2; i++)
    res += v[i];
  return res;
}

inline double prod(const Vec &v)
{
  double res = v[0];
  for (int i = 1; i < 2; i++)
    res *= v[i];
  return res;
}

// element-wise min
inline double min_of(const Vec &v)
{
  double res = v[0];
  for (int i = 1; i < 2; i++)
    if (v[i] < res)
      res = v[i];
  return res;
}

// element-wise max
inline double max_of(const Vec &v)
{
  return -min_of(-v);
}

/*
inline Vec reduce(const Vec &v, int D = 1)
{
  Vec res;
  for (int i = 0; i < D; i++)
    res[i] = v[i];
  for (; i < 1; i++)
    res[i] = v[i + 1];
}
 */
/* inline Vec enlarge(const Vec &v, const double &val, int D = 2)
{
  Vec res;
  for (int i = 0; i < D; i++)
    res[i] = v[i];
  res[D] = val;
  for (int i = D; i < 2; i++)
    res[i + 1] = v[i];
  return res;
} */

//====================================================
// 2d and 3d operators

inline Vec clockwise(const Vec &v)
{
  return Vec{v[1], -v[0]};
}

inline Vec counterclockwise(const Vec &v)
{
  return Vec{-v[1], v[0]};
}

inline double cross(const Vec &lhs, const Vec &rhs)
{
  return lhs[0] * rhs[1] - lhs[1] * rhs[0];
}
