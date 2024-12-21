#pragma once

#include "Vec.h"
#include "Grid.h"
#include "MultiIndexSet.h"
#include "../lib/Eigen/Dense"
#include "../lib/Eigen/Sparse"

using VecList = std::vector<Vec>;
using DiffMatrix = Eigen::SparseMatrix<double>;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

/// @brief Proxy class for Vec, used to set value for VectorValue from a Vec, for there is no Vec &operator() in VectorValue
class VecProxy;

/// @brief Cell average data
class VectorData;
class ScalarData;
/// @brief Boundary face average data
class ScalarBoundaryFaceAvr;
class VectorBoundaryFaceAvr;

/// @brief Face average data
class ScalarFaceAvrData;
class VectorFaceAvrData;









class VecProxy 
{
public:
    VecProxy(double &x, double &y) : x(x), y(y) {}

    VecProxy &operator=(const Vec &vec) 
    {
        x = vec[0];
        y = vec[1];
        return *this;
    }

    operator Vec() const {return Vec{x, y};}

private:
    double &x;
    double &y;
};

class ScalarData
{
public:
    ScalarData() = default;
    ScalarData(const Grid &grid) : grid(grid) 
    {
        data = VectorXd::Zero(grid.get_size()[0] * grid.get_size()[1]);
    }

    ScalarData(const Grid &grid, const VectorXd &data) : grid(grid), data(data) {}
    ScalarData(const ScalarData &other) : grid(other.grid), data(other.data) {}

    ScalarData &operator=(const ScalarData &other) 
    {
        if (this != &other) 
        {
            this->data = other.data;
        }
        return *this;
    }
    //Set and get data
    void set_data(const VectorXd &data) { this->data = data; }
    void set_zero() { data.setZero(); }
    void set_one() { data.setOnes(); }

    const VectorXd &get_data() const { return data; }
    VectorXd &get_data() { return data; }

    friend std::ostream &operator<<(std::ostream &os, const ScalarData &data)
    {
        for (int j=0; j<data.grid.get_size()[1]; j++)
        {
            for (int i=0; i<data.grid.get_size()[0]; i++)
            {
                os << data(i,j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

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
    //Pointwise product
    ScalarData pointwise_product(const ScalarData &rhs) const { return ScalarData(grid, data.cwiseProduct(rhs.get_data())); }
    void pointwise_product(const ScalarData &rhs) { data = data.cwiseProduct(rhs.get_data()); }

    ~ScalarData() = default;
protected:
    const Grid &grid;
    VectorXd data;    
};
ScalarData Multiply(const DiffMatrix &A, const ScalarData &x)
{
    ScalarData result = x;
    result.get_data() = A * x.get_data();
    return result;
}



class VectorData
{
public:
    //Constructors
    VectorData(const Grid &grid) : grid(grid), data_x(grid), data_y(grid) {}
    VectorData(const Grid &grid, const ScalarData &data_x, const ScalarData &data_y) : grid(grid), data_x(data_x), data_y(data_y) {}
    VectorData(const Grid &grid, const VectorXd &data_x, const VectorXd &data_y) : grid(grid), data_x(grid, data_x), data_y(grid, data_y) {}
    VectorData(const VectorData &other) : grid(other.grid), data_x(other.data_x), data_y(other.data_y) {}
    //Set and get data
    VectorData &operator=(const VectorData &other) 
    {
        if (this != &other) 
        {
            this->data_x = other.data_x;
            this->data_y = other.data_y;
        }
        return *this;
    }
    void set_data(const VectorXd &data_x, const VectorXd &data_y) { this->data_x.set_data(data_x); this->data_y.set_data(data_y); }
    const ScalarData &get_data_x() const { return data_x; }
    const ScalarData &get_data_y() const { return data_y; }

    ScalarData &get_data_x() { return data_x; }
    ScalarData &get_data_y() { return data_y; }

    void set_value(int i, int j, const Vec &vec) { data_x(i,j) = vec[0]; data_y(i,j) = vec[1]; }
    void set_value(MultiIndex index, const Vec &vec) { data_x(index) = vec[0]; data_y(index) = vec[1]; }

    void set_zero() { data_x.set_zero(); data_y.set_zero(); }
    void set_one() { data_x.set_one(); data_y.set_one(); }

    Vec operator()(int i, int j) const { return Vec{data_x(i,j), data_y(i,j)}; }
    Vec operator()(MultiIndex index) const { return Vec{data_x(index), data_y(index)}; }

    //Proxy to set value
    VecProxy operator()(int i, int j) { return VecProxy{data_x(i,j), data_y(i,j)}; }
    VecProxy operator()(MultiIndex index) { return VecProxy{data_x(index), data_y(index)}; }
    //Binary operators
    VectorData operator+(const VectorData &rhs) const { return VectorData(grid, data_x + rhs.get_data_x(), data_y + rhs.get_data_y()); }
    VectorData operator-(const VectorData &rhs) const { return VectorData(grid, data_x - rhs.get_data_x(), data_y - rhs.get_data_y()); }
    VectorData operator*(double scalar) const { return VectorData(grid, data_x * scalar, data_y * scalar); }
    VectorData operator/(double scalar) const { return VectorData(grid, data_x / scalar, data_y / scalar); }
    VectorData operator*(Vec vec) const { return VectorData(grid, data_x * vec[0], data_y * vec[1]); }
    VectorData operator/(Vec vec) const { return VectorData(grid, data_x / vec[0], data_y / vec[1]); }

    VectorData Multiply(const DiffMatrix &G_x, const DiffMatrix &G_y) const
    {
        VectorData result(grid);
        result.set_data(G_x*data_x.get_data(), G_y*data_y.get_data());
        return result;
    }


    //Unary operators   
    VectorData &operator+=(const VectorData &rhs) { data_x += rhs.get_data_x(); data_y += rhs.get_data_y(); return *this; }
    VectorData &operator-=(const VectorData &rhs) { data_x -= rhs.get_data_x(); data_y -= rhs.get_data_y(); return *this; }
    VectorData &operator*=(double scalar) { data_x *= scalar; data_y *= scalar; return *this; }
    VectorData &operator/=(double scalar) { data_x /= scalar; data_y /= scalar; return *this; }
    VectorData &operator*=(Vec vec) { data_x *= vec[0]; data_y *= vec[1]; return *this; }
    VectorData &operator/=(Vec vec) { data_x /= vec[0]; data_y /= vec[1]; return *this; }
    //Return the size of the vector data
    int size() const { return data_x.size(); }
    //Return the size of the grid
    int grid_size() const { return grid.get_size()[0] * grid.get_size()[1]; }

    friend std::ostream &operator<<(std::ostream &os, const VectorData &data)
    {
        for (int j=0; j<data.grid.get_size()[1]; j++)
        {
            for (int i=0; i<data.grid.get_size()[0]; i++)
            {
                os << data(i,j) << " ";
            }
            os << std::endl;
        }
        return os;
    }


protected:
    const Grid &grid;
    ScalarData data_x;
    ScalarData data_y;
};

// class VectorData
// {
// public:
//     VectorData(const Grid &grid) : grid(grid) 
//     {
//         data_x = VectorXd::Zero(grid.get_size()[0] * grid.get_size()[1]);
//         data_y = VectorXd::Zero(grid.get_size()[0] * grid.get_size()[1]);
//     }
//     VectorData(const Grid &grid, const VectorXd &data_x, const VectorXd &data_y) : grid(grid), data_x(data_x), data_y(data_y) {}
//     //Set and get data
//     void set_data(const VectorXd &data_x, const VectorXd &data_y) { this->data_x = data_x; this->data_y = data_y; }
//     const VectorXd &get_data_x() const { return data_x; }
//     const VectorXd &get_data_y() const { return data_y; }

//     VectorXd &get_data_x() { return data_x; }
//     VectorXd &get_data_y() { return data_y; }

//     void set_zero() { data_x.setZero(); data_y.setZero(); }
//     void set_one() { data_x.setOnes(); data_y.setOnes(); }
//     //(i,j) or MultiIndex
//     VecProxy operator()(int i, int j) 
//     {
//         int index = grid.MultiToSingle(i, j);
//         return VecProxy{data_x[index], data_y[index]};
//     }
//     VecProxy operator()(MultiIndex index) 
//     {
//         int single_index = grid.MultiToSingle(index);
//         return VecProxy{data_x[single_index], data_y[single_index]};
//     }
//     // Vec &operator()(MultiIndex index) { return Vec{data_x(grid.MultiToSingle(index)), data_y(grid.MultiToSingle(index))}; }

//     //Binary operators
//     VectorData operator+(const VectorData &rhs) const { return VectorData(grid, data_x + rhs.get_data_x(), data_y + rhs.get_data_y()); }
//     VectorData operator-(const VectorData &rhs) const { return VectorData(grid, data_x - rhs.get_data_x(), data_y - rhs.get_data_y()); }
//     VectorData operator*(double scalar) const { return VectorData(grid, data_x * scalar, data_y * scalar); }
//     VectorData operator/(double scalar) const { return VectorData(grid, data_x / scalar, data_y / scalar); }
//     //Unary operators
//     VectorData &operator+=(const VectorData &rhs) { data_x += rhs.get_data_x(); data_y += rhs.get_data_y(); return *this; }
//     VectorData &operator-=(const VectorData &rhs) { data_x -= rhs.get_data_x(); data_y -= rhs.get_data_y(); return *this; }
//     VectorData &operator*=(double scalar) { data_x *= scalar; data_y *= scalar; return *this; }

//     ~VectorData() = default;
// protected:
//     const Grid &grid;
//     VectorXd data_x;
//     VectorXd data_y;   
// };





class ScalarBoundaryFaceAvr
{
public:
    ScalarBoundaryFaceAvr() = default;
    ScalarBoundaryFaceAvr(const Grid &grid) : grid(grid)
    {
        face_avr_l = VectorXd::Zero(grid.get_size()[1]);
        face_avr_r = VectorXd::Zero(grid.get_size()[1]);
        face_avr_d = VectorXd::Zero(grid.get_size()[0]);
        face_avr_u = VectorXd::Zero(grid.get_size()[0]);
    }
    //copy
    ScalarBoundaryFaceAvr(const ScalarBoundaryFaceAvr &other) : grid(other.grid), face_avr_l(other.face_avr_l), face_avr_r(other.face_avr_r), face_avr_d(other.face_avr_d), face_avr_u(other.face_avr_u) {}

    double operator()(int i, int j, int n1, int n2) const;
    double operator()(MultiIndex index, Normal normal) const { return operator()(index[0], index[1], normal[0], normal[1]); }

    double& operator()(int i, int j, int n1, int n2);
    double& operator()(MultiIndex index, Normal normal) { return operator()(index[0], index[1], normal[0], normal[1]); }

    void set_zero() { face_avr_l.setZero(); face_avr_r.setZero(); face_avr_d.setZero(); face_avr_u.setZero(); }
    void set_one() { face_avr_l.setOnes(); face_avr_r.setOnes(); face_avr_d.setOnes(); face_avr_u.setOnes(); }

    VectorXd &get_face_avr_l() { return face_avr_l; }
    VectorXd &get_face_avr_r() { return face_avr_r; }
    VectorXd &get_face_avr_d() { return face_avr_d; }
    VectorXd &get_face_avr_u() { return face_avr_u; }

    const VectorXd &get_face_avr_l() const { return face_avr_l; }
    const VectorXd &get_face_avr_r() const { return face_avr_r; }
    const VectorXd &get_face_avr_d() const { return face_avr_d; }
    const VectorXd &get_face_avr_u() const { return face_avr_u; }

    ScalarBoundaryFaceAvr &operator=(const ScalarBoundaryFaceAvr &other)
    {
        if (this != &other)
        {
            this->face_avr_l = other.face_avr_l;
            this->face_avr_r = other.face_avr_r;
            this->face_avr_d = other.face_avr_d;
            this->face_avr_u = other.face_avr_u;
        }
        return *this;
    }
protected:
    const Grid &grid;
    VectorXd face_avr_l;
    VectorXd face_avr_r;
    VectorXd face_avr_d;
    VectorXd face_avr_u;
};


class VectorBoundaryFaceAvr
{
public:
    VectorBoundaryFaceAvr(const Grid &grid): grid(grid), face_avr_x(grid), face_avr_y(grid) {}
    VectorBoundaryFaceAvr(const Grid &grid, const ScalarBoundaryFaceAvr &face_avr_x, const ScalarBoundaryFaceAvr &face_avr_y) : grid(grid), face_avr_x(face_avr_x), face_avr_y(face_avr_y) {}
    //copy
    VectorBoundaryFaceAvr(const VectorBoundaryFaceAvr &other) : grid(other.grid), face_avr_x(other.face_avr_x), face_avr_y(other.face_avr_y) {}

    void set_data(const ScalarBoundaryFaceAvr &face_avr_x, const ScalarBoundaryFaceAvr &face_avr_y) { this->face_avr_x = face_avr_x; this->face_avr_y = face_avr_y; }
    const ScalarBoundaryFaceAvr &get_face_avr_x() const { return face_avr_x; }
    const ScalarBoundaryFaceAvr &get_face_avr_y() const { return face_avr_y; }

    ScalarBoundaryFaceAvr &get_face_avr_x() { return face_avr_x; }
    ScalarBoundaryFaceAvr &get_face_avr_y() { return face_avr_y; }

    Vec operator()(int i, int j, int n1, int n2) const { return Vec{face_avr_x(i,j,n1,n2), face_avr_y(i,j,n1,n2)}; }
    Vec operator()(MultiIndex index, Normal normal) const { return Vec{face_avr_x(index, normal), face_avr_y(index, normal)}; }

    VecProxy operator()(int i, int j, int n1, int n2) { return VecProxy{face_avr_x(i,j,n1,n2), face_avr_y(i,j,n1,n2)}; }
    VecProxy operator()(MultiIndex index, Normal normal) { return VecProxy{face_avr_x(index, normal), face_avr_y(index, normal)}; }

    VectorBoundaryFaceAvr &operator=(const VectorBoundaryFaceAvr &other)
    {
        if (this != &other)
        {
            this->face_avr_x = other.face_avr_x;
            this->face_avr_y = other.face_avr_y;
        }
        return *this;
    }
protected:
    const Grid &grid;
    ScalarBoundaryFaceAvr face_avr_x;
    ScalarBoundaryFaceAvr face_avr_y;
};



class ScalarFaceAvrData
{
public:
    ScalarFaceAvrData(const Grid &grid) : grid(grid) 
    {
        //Every column denotes
        data_x = MatrixXd::Zero(grid.get_size()[0]+2, grid.get_size()[1]+1);
        data_y = MatrixXd::Zero(grid.get_size()[0]+1, grid.get_size()[1]+2);
    }

    ScalarFaceAvrData(const Grid &grid, const MatrixXd &data_x, const MatrixXd &data_y) : grid(grid), data_x(data_x), data_y(data_y) {}
    //copy
    ScalarFaceAvrData(const ScalarFaceAvrData &other) : grid(other.grid), data_x(other.data_x), data_y(other.data_y) {}


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

    ScalarFaceAvrData &operator=(const ScalarFaceAvrData &other)
    {
        if (this != &other)
        {
            this->data_x = other.data_x;
            this->data_y = other.data_y;
        }
        return *this;
    }

protected:
    const Grid &grid;
    MatrixXd data_x;
    MatrixXd data_y;
};


class VectorFaceAvrData
{
public:
    VectorFaceAvrData(const Grid &grid) : grid(grid), data_x(grid), data_y(grid) {}
    VectorFaceAvrData(const Grid &grid, const ScalarFaceAvrData &data_x, const ScalarFaceAvrData &data_y) : grid(grid), data_x(data_x), data_y(data_y) {}
    //copy
    VectorFaceAvrData(const VectorFaceAvrData &other) : grid(other.grid), data_x(other.data_x), data_y(other.data_y) {}

    void set_data(const ScalarFaceAvrData &data_x, const ScalarFaceAvrData &data_y) { this->data_x = data_x; this->data_y = data_y; }
    const ScalarFaceAvrData &get_data_x() const { return data_x; }
    const ScalarFaceAvrData &get_data_y() const { return data_y; }

    ScalarFaceAvrData &get_data_x() { return data_x; }
    ScalarFaceAvrData &get_data_y() { return data_y; }

    Vec operator()(int i, int j, int n1, int n2) const { return Vec{data_x(i,j,n1,n2), data_y(i,j,n1,n2)}; }
    Vec operator()(MultiIndex index, Normal normal) const { return Vec{data_x(index, normal), data_y(index, normal)}; }

    VecProxy operator()(int i, int j, int n1, int n2) { return VecProxy{data_x(i,j,n1,n2), data_y(i,j,n1,n2)}; }
    VecProxy operator()(MultiIndex index, Normal normal) { return VecProxy{data_x(index, normal), data_y(index, normal)}; }

    VectorFaceAvrData &operator=(const VectorFaceAvrData &other)
    {
        if (this != &other)
        {
            this->data_x = other.data_x;
            this->data_y = other.data_y;
        }
        return *this;
    }
protected:
    const Grid &grid;
    ScalarFaceAvrData data_x;
    ScalarFaceAvrData data_y;
};


double ScalarFaceAvrData::operator()(int i, int j, int n1, int n2) const
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

double& ScalarFaceAvrData::operator()(int i, int j, int n1, int n2)
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


double ScalarBoundaryFaceAvr::operator()(int i, int j, int n1, int n2) const
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

double& ScalarBoundaryFaceAvr::operator()(int i, int j, int n1, int n2)
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