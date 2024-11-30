#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <assert.h>
#include "Vec.h"

using MultiIndex = std::array<int, 2>;

class Grid
{
public:
  // constructors
  Grid();
  Grid(const Vec &_lo, const Vec &_hi, double _h);
  Grid(const Grid &rhs);
  Grid &operator=(const Grid &rhs);

  // accessors
  Vec lo() const;
  Vec hi() const;
  double get_h() const;
  const int *get_size() const;

  // return the index of cell, counting from left-down, starting from 0
  int MultiToSingle(int i, int j) const;
  MultiIndex SingleToMulti(int index) const;


  Vec operator()(int i, int j) const;

  // Based on the left-down coord(i,j)
  Vec center(int i, int j) const;

  friend std::ostream &operator<<(std::ostream &os, const Grid &g);

public:
  bool empty() const;
  bool contain(const Vec &pos) const;
  bool contain(const Grid &rhs) const;
  int volume() const;
  MultiIndex LocateCell(const Vec &pos) const;
  Grid refine() const;
  Grid coarsen() const;

protected:
  Vec corner[2];
  double h;
  int size[2];
};

#define loop_grid_2(gd, i0, i1)                  \
  for (int i1 = 0; i1 <= gd.get_size()[1]; i1++) \
    for (int i0 = 0; i0 <= gd.get_size()[0]; i0++)

#define loop_cell_2(gd, i0, i1)                 \
  for (int i1 = 0; i1 < gd.get_size()[1]; i1++) \
    for (int i0 = 0; i0 < gd.get_size()[0]; i0++)
Grid::Grid() = default;

Grid::Grid(const Vec &_lo, const Vec &_hi, double _h) : corner{_lo, _hi}, h(_h)
{
  for (int i = 0; i < 2; i++)
    size[i] = static_cast<int>(std::ceil((_hi[i] - _lo[i]) / _h));
}

Grid::Grid(const Grid &rhs)
{
  corner[0] = rhs.lo();
  corner[1] = rhs.hi();
  h = rhs.get_h();
  for (int i = 0; i < 2; ++i)
  {
    size[i] = rhs.get_size()[i];
  }
}

Grid &Grid::operator=(const Grid &rhs)
{
  corner[0] = rhs.lo();
  corner[1] = rhs.hi();
  h = rhs.get_h();
  for (int i = 0; i < 2; ++i)
  {
    size[i] = rhs.get_size()[i];
  }
  return *this;
}

Vec Grid::lo() const { return corner[0]; }

Vec Grid::hi() const { return corner[1]; }

double Grid::get_h() const { return h; }

const int *Grid::get_size() const { return size; }

int Grid::MultiToSingle(int i, int j) const
{
  return i + j * size[1];
}

MultiIndex Grid::SingleToMulti(int index) const
{
    int i, j;
    j = index / size[1];
    i = index % size[1];
    return {i, j};
}

Vec Grid::operator()(int i, int j) const
{
  assert(i <= size[0] && j <= size[1]);
  return corner[0] + Vec{i * h, j * h};
}

Vec Grid::center(int i, int j) const
{
  assert(i < size[0] && j < size[1]);
  return corner[0] + Vec{(i + 1. / 2) * h, (j + 1. / 2) * h};
}

std::ostream &operator<<(std::ostream &os, const Grid &g)
{
  for (int i1 = g.size[1]; i1 >= 0; i1--)
  {
    for (int i0 = 0; i0 <= g.size[0]; i0++)
    {
      os << g(i0, i1) << " ";
    }
    os << std::endl;
  }
  return os;
}

bool Grid::empty() const
{
  return min_of(hi() - lo()) < 0;
}

bool Grid::contain(const Vec &pos) const
{
  return (min_of(hi() - pos) >= 0) && (min_of(pos - lo()) >= 0);
}

bool Grid::contain(const Grid &rhs) const
{
  return contain(rhs.lo()) && contain(rhs.hi());
}

int Grid::volume() const
{
  return prod((hi() - lo()) / h);
}

MultiIndex Grid::LocateCell(const Vec &pos) const
{
  MultiIndex index;
  for (int i = 0; i < 2; i++)
  {
    index[i] = static_cast<int>((pos[i] - lo()[i]) / h);
  }
  return index;
}

Grid Grid::refine() const { return Grid(lo(), hi(), h / 2); }

Grid Grid::coarsen() const { return Grid(lo(), hi(), h * 2); }





/** 
class Grid
{
public:
  // constructors
  Grid() = default;

  Grid(const Vec &_lo, const Vec &_hi, double _h) : corner{_lo, _hi}, h(_h)
  {
    for (int i = 0; i < 2; i++)
      size[i] = static_cast<int>(std::ceil((_hi[i] - _lo[i]) / _h));
  }

  Grid(const Grid &rhs)
  {
    corner[0] = rhs.lo();
    corner[1] = rhs.hi();
    h = rhs.get_h();
    for (int i = 0; i < 2; ++i)
    {
      size[i] = rhs.get_size()[i];
    }
  }

  Grid &operator=(const Grid &rhs)
  {
    corner[0] = rhs.lo();
    corner[1] = rhs.hi();
    h = rhs.get_h();
    for (int i = 0; i < 2; ++i)
    {
      size[i] = rhs.get_size()[i];
    }
    return *this;
  }

  // accessors

  Vec lo() const { return corner[0]; }

  Vec hi() const { return corner[1]; }

  double get_h() const { return h; }

  const int *get_size() const { return size; }

  //return the index of cell, couting from left-down, starting from 0
  
  int CellIndex(int i, int j) const
  {
    return i + j * size[1];
  }
  //Verify if the cell is on the boundary

  bool IfOnLeftBound(int i, int j) const
  {
    return i == 0;
  }
  bool IfOnRightBound(int i, int j) const
  {
    return i == size[0] -1;
  }
  bool IfOnDownBound(int i, int j) const
  {
    return j == 0;
  }
  bool IfOnUpBound(int i, int j) const
  {
    return j == size[1] -1;
  }
  //Verify if the cell is near the boundary
  bool IfNearLeftBound(int i, int j) const
  {
    return i == 1;
  }
  bool IfNearRightBound(int i, int j) const
  {
    return i == size[0] -2;
  }
  bool IfNearDownBound(int i, int j) const
  {
    return j == 1;
  }
  bool IfNearUpBound(int i, int j) const
  {
    return j == size[1] -2;
  }

  Vec operator()(int i, int j) const
  {
    assert(i <= size[0] && j <= size[1]);
    return corner[0] + Vec{i * h, j * h};
  }

  // Based on the left-down coord(i,j)
  Vec center(int i, int j) const
  {
    assert(i < size[0] && j < size[1]);
    return corner[0] + Vec{(i + 1. / 2) * h, (j + 1. / 2) * h};
  }

  friend std::ostream &operator<<(std::ostream &os, const Grid &g)
  {
    for (int i1 = g.size[1]; i1 >= 0; i1--)
    {
      for (int i0 = 0; i0 <= g.size[0]; i0++)
      {
        os << g(i0, i1) << " ";
      }
      os << std::endl;
    }
    return os;
  }

public:
  bool empty() const
  {
    return min_of(hi() - lo()) < 0;
  }

  bool contain(const Vec &pos) const
  {
    return (min_of(hi() - pos) >= 0) && (min_of(pos - lo()) >= 0);
  }

  bool contain(const Grid &rhs) const
  {
    return contain(rhs.lo()) && contain(rhs.hi());
  }

  int volume() const
  {
    return prod((hi() - lo()) / h);
  }

  Grid refine() const { return Grid(lo(), hi(), h / 2); }

  Grid coarsen() const { return Grid(lo(), hi(), h * 2); }

protected:
  Vec corner[2];
  double h;
  int size[2];
};

#define loop_grid_2(gd, i0, i1)                  \
  for (int i1 = 0; i1 <= gd.get_size()[1]; i1++) \
    for (int i0 = 0; i0 <= gd.get_size()[0]; i0++)

#define loop_cell_2(gd, i0, i1)                 \
  for (int i1 = 0; i1 < gd.get_size()[1]; i1++) \
    for (int i0 = 0; i0 < gd.get_size()[0]; i0++)
    **/