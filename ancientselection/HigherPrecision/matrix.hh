#ifndef MATRIX_HH
#define MATRIX_HH

#include <fstream>
#include "matrix.h"

template<typename T, typename U>
void Matrix<T,U>::allocate() {
    rawData = new T [nx*ny];
    grid    = new T* [nx];
    for (U iX=0; iX<nx; ++iX) {
        grid[iX] = rawData + iX*ny;
    }
}

template<typename T, typename U>
Matrix<T,U>::Matrix(U nx_, U ny_) : nx(nx_), ny(ny_) {
    allocate();
}

template<typename T, typename U>
Matrix<T,U>::Matrix(const Matrix<T,U> &rhs) : nx(rhs.getNx()), ny(rhs.getNy()) {
    allocate();
    for (U iX = 0; iX < nx; ++iX) {
        for (U iY = 0; iY < ny; ++iY) {
            grid[iX][iY] = rhs.get(iX,iY);
        }
    }
}

template<typename T, typename U>
Matrix<T,U>::Matrix(U N) : nx(N), ny(N) {
    allocate();
}

template<typename T, typename U>
Matrix<T,U>::~Matrix() {
    delete[] rawData;
    delete[] grid;
}

template<typename T, typename U>
Matrix<T,U>& Matrix<T,U>::operator= (Matrix<T,U> const& rhs ) {
    Matrix<T,U> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, typename U>
void Matrix<T,U>::swap(Matrix<T,U> &rhs) {
    std::swap(rawData, rhs.rawData);
    std::swap(grid, rhs.grid);
}

template<typename T, typename U>
void Matrix<T,U>::print(int len) const {
    for (U iX = 0; iX < getNx(); ++iX) {
        for (U iY = 0; iY < getNy(); ++iY) {
            std::cout << std::setprecision(len) << get(iX,iY) << " ";
        }
        std::cout << std::setprecision(len) <<  std::endl;
    }
}

template<typename T, typename U>
void Matrix<T,U>::print(std::string fname, int len) const {
    std::ofstream fout(fname.c_str());
    
    for (U iX = 0; iX < getNx(); ++iX) {
        for (U iY = 0; iY < getNy(); ++iY) {
            fout << std::setprecision(len) << get(iX,iY) << " ";
        }
        fout << std::setprecision(len) <<  std::endl;
    }
}

#endif // MATRIX_HH
