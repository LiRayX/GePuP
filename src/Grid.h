#pragma once

#include <iostream>
#include <assert.h>
#include <array>
#include "Vec.h"


using MultiIndex = std::array<int, 2>;
using Normal = std::array<int, 2>;
using VecList = std::vector<Vec>;
using MultiIndexList = std::vector<MultiIndex>;


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

  double get_cell_volume() const {return h*h;}
  // return the index of cell, counting from left-down, starting from 0
  int MultiToSingle(int i, int j) const;
  int MultiToSingle(MultiIndex index) const;


  MultiIndex SingleToMulti(int index) const;

  MultiIndexList getOuterCells() const;

  int layer(const MultiIndex &ghost_index) const;
  //Check if the Point is on the Face
  bool OnFace(int i, int j, int x, int y, const Vec &point, double tol) const;
  bool OnFace(MultiIndex index, Normal normal, const Vec &point, double tol) const;

  //Get 4 corner points of the cell
  VecList getAllCorners(const MultiIndex &index) const;

  //Check if the Point is on the Corner
  bool isCorner(const MultiIndex &index, const Vec &point, double tol) const;
  //Given two intersection points on the vertical faces, return the corner point
  Vec getCutCorner(const MultiIndex &index, const Vec &intersection_1, const Vec &intersection_2, double tol) const;
  //Given two intersection points on the parallel face, return the corner points
  
  //Sign distance from the point to the face
  double SignDistance(MultiIndex index, Normal normal,const Vec &point) const;
  double Distance(MultiIndex index, Normal normal,const Vec &point) const;
  //Specify the face of the point located
  Normal getDirection(MultiIndex index, const Vec &point, double tol) const;

  //corner points located in the left-down side and the index should valids
  Vec operator()(int i, int j) const;
  Vec operator()(MultiIndex index) const;
  //Return the center of the cell
  // Based on the left-down coord(i,j)
  Vec center(int i, int j) const;
  Vec center(MultiIndex index) const;

  //Check if the index is valid
  bool isIndexValid(int i, int j) const;
  bool isIndexValid(MultiIndex index) const;


  friend std::ostream &operator<<(std::ostream &os, const Grid &g);

public:
  bool empty() const;
  bool contain(const Vec &pos) const;
  bool contain(const Grid &rhs) const;
  // return the volume of the whole domain
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

#define loop_inner_cell_2(gd, i0, i1)             \
  for (int i1 = 2; i1 < gd.get_size()[1]-2; i1++) \
    for (int i0 = 2; i0 < gd.get_size()[0]-2; i0++)

                              

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
int Grid::MultiToSingle(MultiIndex index) const
{
  return index[0] + index[1] * size[1];
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

Vec Grid::operator()(MultiIndex index) const
{
  return operator()(index[0], index[1]);
}

Vec Grid::center(int i, int j) const
{
  assert(i < size[0] && j < size[1]);
  return corner[0] + Vec{(i + 1. / 2) * h, (j + 1. / 2) * h};
}

Vec Grid::center(MultiIndex index) const
{
  return center(index[0], index[1]);
}

VecList Grid::getAllCorners(const MultiIndex &index) const
{
  VecList corners;
  corners.push_back(center(index)+Vec{-1/2.0*h, -1/2.0*h});
  corners.push_back(center(index)+Vec{1/2.0*h, -1/2.0*h});
  corners.push_back(center(index)+Vec{1/2.0*h, 1/2.0*h});
  corners.push_back(center(index)+Vec{-1/2.0*h, 1/2.0*h});
  return corners;
}

bool Grid::OnFace(int i, int j, int x, int y, const Vec &point, double tol) const
{
  assert((x == 0 && (y == 1 || y == -1)) || (y == 0 && (x == 1 || x == -1)));
  Vec cell_center = this->center(i, j);
  if(x==0)
  {
    return (std::fabs(point[1]-(cell_center[1] + y/2.0 *h)) < tol)&&(std::fabs(point[0]-cell_center[0]) < 1/2.0 *h);
  }
  else
  {
    return (std::fabs(point[0]-(cell_center[0] + x/2.0 *h)) < tol)&&(std::fabs(point[1]-cell_center[1]) < 1/2.0 *h);
  }
}
bool Grid::OnFace(MultiIndex index, Normal normal, const Vec &point, double tol) const
{
  assert((normal[0] == 0 && (normal[1] == 1 || normal[1] == -1)) || (normal[1] == 0 && (normal[0] == 1 || normal[0] == -1)));
  Vec cell_center = this->center(index[0], index[1]);
  if(normal[0]==0)
  {
    return (std::fabs(point[1]-(cell_center[1] + normal[1]/2.0 *h)) < tol)&&(std::fabs(point[0]-cell_center[0]) < 1/2.0 *h);
  }
  else
  {
    return (std::fabs(point[0]-(cell_center[0] + normal[0]/2.0 *h)) < tol)&&(std::fabs(point[1]-cell_center[1]) < 1/2.0 *h);
  }
}
double Grid::SignDistance(MultiIndex index, Normal normal, const Vec &point) const
{
  assert((normal[0] == 0 && (normal[1] == 1 || normal[1] == -1)) || (normal[1] == 0 && (normal[0] == 1 || normal[0] == -1)));
  Vec cell_center = this->center(index[0], index[1]);
  if(normal[0] == 0)
  {
    assert(std::fabs(point[0] - cell_center[0]) < 3/2.0 * h);
    double diff = point[1] - (cell_center[1] + normal[1] / 2.0 * h);
    // return (std::fabs(diff) >= tol) ? diff : 0;
    return diff;
  }
  else
  {
    assert(std::fabs(point[1] - cell_center[1]) < 3/2.0 * h);
    double diff = point[0] - (cell_center[0] + normal[0] / 2.0 * h);
    // return (std::fabs(diff) >= tol) ? diff : 0;
    return diff;
  }
}

double Grid::Distance(MultiIndex index, Normal normal, const Vec &point) const
{
  return std::fabs(SignDistance(index, normal, point));
}

Normal Grid::getDirection(MultiIndex index, const Vec &point, double tol) const
{
  Normal normal;
  if (fabs(fabs(point[0] - center(index)[0]) - 1.0/2*h) < tol)
  {
    int sign = (point[0] - center(index)[0] > 0) ? 1 : -1;
    return Normal{sign, 0};
  }
  else if (fabs(fabs(point[1] - center(index)[1]) - 1.0/2*h) < tol)
  {
    int sign = (point[1] - center(index)[1] > 0) ? 1 : -1;
    return Normal{0, sign};
  }
  else
  {
    assert(false);
  }

  // if (Distance(index, {0, 1}, point) < tol)
  // {
  //   normal = {0, 1};
  // }
  // else if (Distance(index, {0, -1}, point) < tol)
  // {
  //   normal = {0, -1};
  // }
  // else if (Distance(index, {1, 0}, point) < tol)
  // {
  //   normal = {1, 0};
  // }
  // else if (Distance(index, {-1, 0}, point) < tol)
  // {
  //   normal = {-1, 0};
  // }
  // else
  // {
  //   normal = {0, 0};
  // }
}


bool Grid::isCorner(const MultiIndex &index, const Vec &point, double tol = 1e-12) const
{
  if (norm(abs(point-center(index)) - Vec{0.5*h, 0.5*h}) < tol)
  {
    return true;
  }
  else
  {
    return false;
  }
}

Vec Grid::getCutCorner(const MultiIndex &index, const Vec &intersection_1, const Vec &intersection_2, double tol = 1e-12) const
{
  Vec cut_corner{intersection_1[0], intersection_2[1]};
  if (isCorner(index, cut_corner, tol))
  {
    return cut_corner;
  }
  else if(isCorner(index, {intersection_2[0], intersection_1[1]}, tol))
  {
    cut_corner = {intersection_2[0], intersection_1[1]};
    return cut_corner;
  }
  else
  {
    assert(false);
  }
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

bool Grid::isIndexValid(int i, int j) const
{
  return i >= 0 && i < size[0] && j >= 0 && j < size[1];
}
bool Grid::isIndexValid(MultiIndex index) const
{
  return isIndexValid(index[0], index[1]);
}

MultiIndexList Grid::getOuterCells() const
{
  MultiIndexList outerCells;
  int m = size[0];
  int n = size[1];
  for (int i1 = 0; i1 < 2; i1++)                  
      for (int i0 = 0; i0 < m; i0++)
          outerCells.push_back({i0, i1}); 
  for (int i1 = n - 2; i1 < n; i1++) 
      for (int i0 = 0; i0 < m; i0++) 
          outerCells.push_back({i0, i1});
  for (int i1 = 0; i1 < 2; i1++)                  
      for (int i0 = 2; i0 < n - 3; i0++)
          outerCells.push_back({i0, i1});
  for (int i1 = n - 2; i1 < n; i1++) 
      for (int i0 = 2; i0 < m - 2; i0++)
          outerCells.push_back({i0, i1});
  return outerCells;
}

Grid Grid::refine() const { return Grid(lo(), hi(), h / 2); }

Grid Grid::coarsen() const { return Grid(lo(), hi(), h * 2); }

int Grid::layer(const MultiIndex &ghost_index) const
{
  int i = ghost_index[0];
  int j = ghost_index[1];
  if (i==-1 || i==size[0] || j==-1 || j==size[1])
  {
    return 1;
  }
  else if (i==-2 || i==size[0]+1 || j==-2 || j==size[1]+1)
  {
    return 2;
  }
  else
  {
    return 0;
  }
}




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