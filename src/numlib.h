#pragma once

inline constexpr int sign(double a)
{
  if (a > 0)
    return 1;
  if (a < 0)
    return -1;
  return 0;
}

inline constexpr int factorial(int n)
{
  int r = 1;
  for (; n > 1; r *= n, --n)
    ;
  return r;
}

inline constexpr int binom(int n, int m)
{
  return factorial(n) / factorial(n - m) / factorial(m);
}

struct GaussLegendreConstant4order
{
  static constexpr double knots[] = {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
  static constexpr double weights[] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
};

template <class T_Func>
inline double quad(const T_Func &g, double a, double b)
{
  double R = 0;
  double u = (a + b) / 2, v = (b - a) / 2;
  int k = 0;
  for (; k < 4; ++k)
    R += GaussLegendreConstant4order::weights[k] * g(u + GaussLegendreConstant4order::knots[k] * v);
  return v * R;
}