#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <cmath>
#include <memory>


template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > prod(Matrix<T,U> &A, Matrix<T,U> &B);

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > prod(const Matrix<T,U> &A, T alpha);

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > sub(const Matrix<T,U> &A, const Matrix<T,U> &B);

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > add(const Matrix<T,U> &A, const Matrix<T,U> &B);

template<typename T, typename U>
T norm(const Matrix<T,U> &A);

template<typename T, typename U>
T trace(const Matrix<T,U> &A);

template<typename T, typename U>
bool ensurseCompatibilityWithProbabilities(Matrix<T,U> &A, T threshold);

template<typename T, typename U>
bool ensursePositivityAndOneSmallness(Matrix<T,U> &A, T threshold);

#endif // MATRIX_OPERATIONS_H
