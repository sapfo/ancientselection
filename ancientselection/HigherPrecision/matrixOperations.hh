#ifndef MATRIX_OPERATIONS_HH
#define MATRIX_OPERATIONS_HH

#include <cmath>

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > prod(Matrix<T,U> &A, Matrix<T,U> &B) {
    assert(A.getNy() == B.getNx());
    Matrix<T,U> *res = new Matrix<T,U>(A.getNx(), B.getNy());
    
//     Rgemm("t", "t", A.getNx(), A.getNx(), A.getNx(), (T)1, A.get(), A.getNx(), B.get(), A.getNx(), T(), res->get(), A.getNx());
//     res->transposeInPlace();
    
    for (U iA = 0; iA < A.getNx(); ++iA) {
        for (U iB = 0; iB < B.getNy(); ++iB) {
            res->get(iA,iB) = T();
            for (U iC = 0; iC < A.getNy(); ++iC) {
                res->get(iA,iB) += A.get(iA,iC)*B.get(iC,iB);
            }
        }
    }
    
    return std::auto_ptr<Matrix<T,U> >(res);
}

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > prod(const Matrix<T,U> &A, const T alpha) {
    Matrix<T,U> *res = new Matrix<T,U>(A.getNx(), A.getNy());
    
    for (U iX = 0; iX < A.getNx(); ++iX) {
        for (U iY = 0; iY < A.getNy(); ++iY) {
            res->get(iX,iY) = alpha * A.get(iX,iY);
        }
    }
    return std::auto_ptr<Matrix<T,U> >(res);
}

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > add(const Matrix<T,U> &A, const Matrix<T,U> &B) {
    assert(A.getNx() == B.getNx());
    assert(A.getNy() == B.getNy());
    
    Matrix<T,U> *res = new Matrix<T,U>(A.getNx(), A.getNy());
    
    for (U iX = 0; iX < A.getNx(); ++iX) {
        for (U iY = 0; iY < A.getNy(); ++iY) {
            res->get(iX,iY) = A.get(iX,iY) + B.get(iX,iY);
        }
    }
    return std::auto_ptr<Matrix<T,U> >(res);
}

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > sub(const Matrix<T,U> &A, const Matrix<T,U> &B) {
    assert(A.getNx() == B.getNx());
    assert(A.getNy() == B.getNy());
    
    Matrix<T,U> *res = new Matrix<T,U>(A.getNx(), A.getNy());
    
    for (U iX = 0; iX < A.getNx(); ++iX) {
        for (U iY = 0; iY < A.getNy(); ++iY) {
            res->get(iX,iY) = A.get(iX,iY) - B.get(iX,iY);
        }
    }
    return std::auto_ptr<Matrix<T,U> >(res);
}

template<typename T, typename U>
T norm(const Matrix<T,U> &A) {
    T tot = T();
    
    for (U iX = 0; iX < A.getNx(); ++iX) {
        for (U iY = 0; iY < A.getNy(); ++iY) {
            tot += abs(A.get(iX,iY));
        }
    }
    return tot;
}

template<typename T, typename U>
T trace(const Matrix<T,U> &A) {
    assert(A.getNx() == A.getNy());
    T tot = T();
    
    for (U iX = 0; iX < A.getNx(); ++iX) {
        tot += A.get(iX,iX);
    }
    return tot;
}

template<typename T, typename U>
bool ensurseCompatibilityWithProbabilities(Matrix<T,U> &A, T threshold) {
    for (U iX = 0; iX < A.getNx(); ++iX) {
        T line = T();
        for (U iY = 0; iY < A.getNy(); ++iY) {
            line += A.get(iX,iY);
            if (A.get(iX,iY) < T()) {
                if (A.get(iX,iY) > -threshold) {
                    A.get(iX,iY) = T();
                }
                else {
                    std::cout << "A smaller than 0." << std::endl;
                    std::cout << A.get(iX,iY) << std::endl;
                    return false;
                }
            }
            else if (A.get(iX,iY) > (T)1) {
                if (A.get(iX,iY)-(T)1 < threshold) {
                    A.get(iX,iY) = (T)1;
                }
                else {
                    std::cout << "A bigger than 1." << std::endl;
                    std::cout << A.get(iX,iY) << std::endl;
                    return false;
                }
            }
        }
        if (line < (T)1-threshold && line > (T)1+threshold) {
            std::cout << "line = " << std::endl;
            std::cout << line << std::endl;
            return false;
        }
    }
    return true;
}

template<typename T, typename U>
bool ensursePositivityAndOneSmallness(Matrix<T,U> &A, T threshold) {
    for (U iX = 0; iX < A.getNx(); ++iX) {
        for (U iY = 0; iY < A.getNy(); ++iY) {
            if (A.get(iX,iY) < T()) {
                if (A.get(iX,iY) > -threshold) {
                    A.get(iX,iY) = T();
                }
                else {
                    std::cout << "A smaller than 0." << std::endl;
                    std::cout << A.get(iX,iY) << std::endl;
                    return false;
                }
            }
            else if (A.get(iX,iY) > (T)1) {
                if (A.get(iX,iY)-(T)1 < threshold) {
                    A.get(iX,iY) = (T)1;
                }
                else {
                    std::cout << "A bigger than 1." << std::endl;
                    std::cout << A.get(iX,iY) << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

#endif // MATRIX_OPERATIONS_HH
