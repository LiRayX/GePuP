#pragma once

#include "Vec.h"
#include "Grid.h"
#include "MultiIndexSet.h"
#include "../lib/Eigen/Dense"
#include "../lib/Eigen/Sparse"

using VecList = std::vector<Vec>;
using SparseMatrix = Eigen::SparseMatrix<double>;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

class ScalarData
{
public:

    ScalarData(const Grid &grid) : grid(grid) 
    {
        data = VectorXd::Zero(grid.get_size()[0] * grid.get_size()[1]);
    }
    ScalarData(const Grid &grid, const VectorXd &data) : grid(grid), data(data) {}
    //Set and get data
    void set_data(const VectorXd &data) { this->data = data; }
    const VectorXd &get_data() const { return data; }
    void set_zero() { data.setZero(); }
    void set_one() { data.setOnes(); }
    //MultiIndex
    double operator()(int i, int j) const { return data(grid.MultiToSingle(i,j)); }
    double operator()(MultiIndex index) const { return data(grid.MultiToSingle(index)); }
    double &operator()(int i, int j) { return data(grid.MultiToSingle(i,j)); }
    double &operator()(MultiIndex index) { return data(grid.MultiToSingle(index)); }
    //Binary operators
    ScalarData operator+(const ScalarData &rhs) const { return ScalarData(grid, data + rhs.get_data()); }
    ScalarData operator-(const ScalarData &rhs) const { return ScalarData(grid, data - rhs.get_data()); }
    ScalarData operator*(double scalar) const { return ScalarData(grid, data * scalar); }
    ScalarData operator/(double scalar) const { return ScalarData(grid, data / scalar); }
    //Unary operators
    ScalarData &operator+=(const ScalarData &rhs) { data += rhs.get_data(); return *this; }
    ScalarData &operator-=(const ScalarData &rhs) { data -= rhs.get_data(); return *this; }
    ScalarData &operator*=(double scalar) { data *= scalar; return *this; }
    ScalarData &operator/=(double scalar) { data /= scalar; return *this; }
    //Return the size of the scalar data
    int size() const { return data.size(); }
    //Return the size of the grid
    int grid_size() const { return grid.get_size()[0] * grid.get_size()[1]; }


    ~ScalarData() = default;
protected:
    const Grid &grid;
    VectorXd data;    
};


class VecProxy {
public:
    VecProxy(double &x, double &y) : x(x), y(y) {}

    VecProxy &operator=(const Vec &vec) {
        x = vec[0];
        y = vec[1];
        return *this;
    }

    operator Vec() const {return Vec{x, y};}

private:
    double &x;
    double &y;
};


class VectorData
{
public:
    VectorData(const Grid &grid) : grid(grid) 
    {
        data_x = VectorXd::Zero(grid.get_size()[0] * grid.get_size()[1]);
        data_y = VectorXd::Zero(grid.get_size()[0] * grid.get_size()[1]);
    }
    VectorData(const Grid &grid, const VectorXd &data_x, const VectorXd &data_y) : grid(grid), data_x(data_x), data_y(data_y) {}
    //Set and get data
    void set_data(const VectorXd &data_x, const VectorXd &data_y) { this->data_x = data_x; this->data_y = data_y; }
    const VectorXd &get_data_x() const { return data_x; }
    const VectorXd &get_data_y() const { return data_y; }

    VectorXd &get_data_x() { return data_x; }
    VectorXd &get_data_y() { return data_y; }

    void set_zero() { data_x.setZero(); data_y.setZero(); }
    void set_one() { data_x.setOnes(); data_y.setOnes(); }
    //(i,j) or MultiIndex
    VecProxy operator()(int i, int j) 
    {
        int index = grid.MultiToSingle(i, j);
        return VecProxy{data_x[index], data_y[index]};
    }
    VecProxy operator()(MultiIndex index) 
    {
        int single_index = grid.MultiToSingle(index);
        return VecProxy{data_x[single_index], data_y[single_index]};
    }
    // Vec &operator()(MultiIndex index) { return Vec{data_x(grid.MultiToSingle(index)), data_y(grid.MultiToSingle(index))}; }

    //Binary operators
    VectorData operator+(const VectorData &rhs) const { return VectorData(grid, data_x + rhs.get_data_x(), data_y + rhs.get_data_y()); }
    VectorData operator-(const VectorData &rhs) const { return VectorData(grid, data_x - rhs.get_data_x(), data_y - rhs.get_data_y()); }
    VectorData operator*(double scalar) const { return VectorData(grid, data_x * scalar, data_y * scalar); }
    VectorData operator/(double scalar) const { return VectorData(grid, data_x / scalar, data_y / scalar); }
    //Unary operators
    VectorData &operator+=(const VectorData &rhs) { data_x += rhs.get_data_x(); data_y += rhs.get_data_y(); return *this; }
    VectorData &operator-=(const VectorData &rhs) { data_x -= rhs.get_data_x(); data_y -= rhs.get_data_y(); return *this; }
    VectorData &operator*=(double scalar) { data_x *= scalar; data_y *= scalar; return *this; }

    ~VectorData() = default;
protected:
    const Grid &grid;
    VectorXd data_x;
    VectorXd data_y;   
};



class BoundaryFaceAvr
{
public:
    BoundaryFaceAvr(const Grid &grid, const MultiIndexSet &CellOnBoundary) : grid(grid), CellOnBoundary(CellOnBoundary) 
    {
        face_avr_l = VectorXd::Zero(grid.get_size()[1]);
        face_avr_r = VectorXd::Zero(grid.get_size()[1]);
        face_avr_d = VectorXd::Zero(grid.get_size()[0]);
        face_avr_u = VectorXd::Zero(grid.get_size()[0]);
    }
    ~BoundaryFaceAvr();

    double operator()(int i, int j, int n1, int n2) const;
    double operator()(MultiIndex index, Normal normal) const { return operator()(index[0], index[1], normal[0], normal[1]); }

    double& operator()(int i, int j, int n1, int n2);
    double& operator()(MultiIndex index, Normal normal) { return operator()(index[0], index[1], normal[0], normal[1]); }

    void set_zero() { face_avr_l.setZero(); face_avr_r.setZero(); face_avr_d.setZero(); face_avr_u.setZero(); }
    void set_one() { face_avr_l.setOnes(); face_avr_r.setOnes(); face_avr_d.setOnes(); face_avr_u.setOnes(); }

    const VectorXd &get_face_avr_l() const { return face_avr_l; }
    const VectorXd &get_face_avr_r() const { return face_avr_r; }
    const VectorXd &get_face_avr_d() const { return face_avr_d; }
    const VectorXd &get_face_avr_u() const { return face_avr_u; }


protected:
    const Grid &grid;
    const MultiIndexSet &CellOnBoundary;
    VectorXd face_avr_l;
    VectorXd face_avr_r;
    VectorXd face_avr_d;
    VectorXd face_avr_u;
};


class FaceAvrData
{
public:
    FaceAvrData(const Grid &grid) : grid(grid) 
    {
        //Every column denotes
        data_x = MatrixXd::Zero(grid.get_size()[0]+2, grid.get_size()[1]+1);
        data_y = MatrixXd::Zero(grid.get_size()[0]+1, grid.get_size()[1]+2);
    }

    FaceAvrData(const Grid &grid, const MatrixXd &data_x, const MatrixXd &data_y) : grid(grid), data_x(data_x), data_y(data_y) {}

    void set_data(const MatrixXd &data_x, const MatrixXd &data_y) { this->data_x = data_x; this->data_y = data_y; }
    const MatrixXd &get_data_x() const { return data_x; }
    const MatrixXd &get_data_y() const { return data_y; }

    MatrixXd &get_data_x() { return data_x; }
    MatrixXd &get_data_y() { return data_y; }

    void set_zero() { data_x.setZero(); data_y.setZero(); }
    void set_one() { data_x.setOnes(); data_y.setOnes(); }

    double operator()(int i, int j, int n1, int n2) const;
    double operator()(MultiIndex index, Normal normal) const { return operator()(index[0], index[1], normal[0], normal[1]); }

    double &operator()(int i, int j, int n1, int n2);
    double &operator()(MultiIndex index, Normal normal) { return operator()(index[0], index[1], normal[0], normal[1]); }


protected:
    const Grid &grid;
    MatrixXd data_x;
    MatrixXd data_y;
};


double FaceAvrData::operator()(int i, int j, int n1, int n2) const
{
    assert((n1 == 0 && (n2 == 1 || n2 == -1)) || (n2 == 0 && (n1 == 1 || n1 == -1)));
    if(n1 == 0 && n2 == 1)
    {
        return data_x(i+1,j+1);
    }
    else if(n1 == 0 && n2 == -1)
    {
        return data_x(i+1,j);
    }
    else if(n1 == 1 && n2 == 0)
    {
        return data_y(i+1,j+1);
    }
    else if(n1 == -1 && n2 == 0)
    {
        return data_y(i,j+1);
    }
    throw std::out_of_range("Invalid indices for FaceAvrData::operator()");
}

double& FaceAvrData::operator()(int i, int j, int n1, int n2)
{
    assert((n1 == 0 && (n2 == 1 || n2 == -1)) || (n2 == 0 && (n1 == 1 || n1 == -1)));
    if(n1 == 0 && n2 == 1)
    {
        return data_x(i+1,j+1);
    }
    else if(n1 == 0 && n2 == -1)
    {
        return data_x(i+1,j);
    }
    else if(n1 == 1 && n2 == 0)
    {
        return data_y(i+1,j+1);
    }
    else if(n1 == -1 && n2 == 0)
    {
        return data_y(i,j+1);
    }
    throw std::out_of_range("Invalid indices for FaceAvrData::operator()");
}


double BoundaryFaceAvr::operator()(int i, int j, int n1, int n2) const
{
    assert((n1 == 0 && (n2 == 1 || n2 == -1)) || (n2 == 0 && (n1 == 1 || n1 == -1)));
    assert(i==0 || i==grid.get_size()[0]-1 || j==0 || j==grid.get_size()[1]-1);
    if(n1 == 0 && n2 == 1)
    {
        return face_avr_u[i];
    }
    else if(n1 == 0 && n2 == -1)
    {
        return face_avr_d[i];
    }
    else if(n1 == 1 && n2 == 0)
    {
        return face_avr_r[j];
    }
    else if(n1 == -1 && n2 == 0)
    {
        return face_avr_l[j];
    }
    throw std::out_of_range("Invalid indices for BoundaryFaceAvr::operator()");
}

double& BoundaryFaceAvr::operator()(int i, int j, int n1, int n2)
{
    assert((n1 == 0 && (n2 == 1 || n2 == -1)) || (n2 == 0 && (n1 == 1 || n1 == -1)));
    assert(i==0 || i==grid.get_size()[0]-1 || j==0 || j==grid.get_size()[1]-1);
    if(n1 == 0 && n2 == 1)
    {
        return face_avr_u[i];
    }
    else if(n1 == 0 && n2 == -1)
    {
        return face_avr_d[i];
    }
    else if(n1 == 1 && n2 == 0)
    {
        return face_avr_r[j];
    }
    else if(n1 == -1 && n2 == 0)
    {
        return face_avr_l[j];
    }
    throw std::out_of_range("Invalid indices for BoundaryFaceAvr::operator()");
}