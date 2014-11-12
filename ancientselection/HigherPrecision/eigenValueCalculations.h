#ifndef EIGEN_VALUE_CALCULATIONS_H
#define EIGEN_VALUE_CALCULATIONS_H

#include "matrix.h"
#include "vector.h"
#include <vector>

template<typename T, typename U>
void eigenValuesSymmetricMatrixUpper(Matrix<T,U> &A, Vector<T,U> &w);

template<typename T, typename U>
void eigenValuesSymmetricMatrixLower(Matrix<T,U> &A, Vector<T,U> &w);

template<typename T, typename U>
void eigenValuesAndVectorsSymmetricMatrix(Matrix<T,U> &A, Vector<T,U> &w);

template<typename T, typename U>
Matrix<T,U> exponentiateMatrix(T time, T gamma, T h, U n);

template<typename T, typename U>
std::vector<Matrix<T,U> > exponentiateBigMatrix(std::vector<T> &timeVec, T gamma, T h, U n);

#endif // EIGENVALUE_CALCULATIONS_H
