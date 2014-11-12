#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

template<typename T, typename U>
class Matrix {
public :
    Matrix() { };
    Matrix(U N);
    Matrix(U nx_, U ny_);
    Matrix(const Matrix<T,U> &rhs);
    ~Matrix();
    
    Matrix<T,U>& operator=(Matrix<T,U> const& rhs);
    void swap(Matrix<T,U> &rhs);
public :
    // x extent of the matrix
    U getNx() const {
        return nx;
    }
    // y extent of the matrix
    U getNy() const {
        return ny;
    }
    // returns the array containing the matrix data
    T *get() {
        return rawData;
    }
    // read-write access to an index of the matrix
    T & get(U iX, U iY) {
        assert(iX < nx);
        assert(iY < ny);
        return grid[iX][iY];
    }
    // read only access to an index of the matrix
    T const& get(U iX, U iY) const {
        assert(iX < nx);
        assert(iY < ny);
        return grid[iX][iY];
    }
    
    // extracts a row of the matrix 
    Vector<T,U> row(U iX) const {
        Vector<T,U> vec(ny);
        for (U iY = 0; iY < ny; ++iY) {
            vec.get(iY) = get(iX,iY);
        }
        return vec;
    }
    
    // extracts a column of the matrix 
    Vector<T,U> col(U iY) const {
        Vector<T,U> vec(nx);
        for (U iX = 0; iX < nx; ++iX) {
            vec.get(iX) = get(iX,iY);
        }
        return vec;
    }
    
    // transposes "in place" the matrix (old values are lost)
    void transposeInPlace() {
        assert(nx == ny);
        for (U iX = 0; iX < nx; ++iX) {
            for (U iY = iX; iY < ny; ++iY) {
                T tmp = get(iX,iY);
                get(iX,iY) = get(iY,iX);
                get(iY,iX) = tmp;
            }
        }
    }
    
    // set all indices of the matrix to the "val" value.
    void setToConst(T val) {
        for (U iX = 0; iX < nx; ++iX) {
            for (U iY = 0; iY < ny; ++iY) {
                get(iX,iY) = val;
            }
        }
    }
    
    void print(int len = 15) const;
    void print(std::string fname, int len = 15) const;
private :
    void allocate();
private :
    U nx, ny;
    T *rawData;
    T **grid;
};

#endif // MATRIX_H
