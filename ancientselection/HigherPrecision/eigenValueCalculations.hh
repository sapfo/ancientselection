#ifndef EIGEN_VALUE_CALCULATIONS_HH
#define EIGEN_VALUE_CALCULATIONS_HH

#include "eigenValueCalculations.h"
#include "constantsHelpers.h"
#include "matrixOperations.h"
#include "matrixOperations.hh"

template<typename T, typename U>
void eigenValuesSymmetricMatrixUpper(Matrix<T,U> &A, Vector<T,U> &w) {
    U n = A.getNx();
    Vector<T,U> workTest(1);
    U lwork = -1, info = 0;
    
    dsyev_( "N", "U", &n, A.get(), &n, w.get(), workTest.get(), &lwork, &info );
    lwork= workTest.get(0);
    Vector<T,U> work(lwork);
    dsyev_( "N", "U", &n, A.get(), &n, w.get(), work.get(), &lwork, &info );
    assert(info == 0);
}

template<typename T, typename U>
void eigenValuesSymmetricMatrixLower(Matrix<T,U> &A, Vector<T,U> &w) {
    U n = A.getNx();
    Vector<T,U> workTest(1);
    U lwork = -1, info = 0;
    
    dsyev_( "N", "L", &n, A.get(), &n, w.get(), workTest.get(), &lwork, &info );
    lwork= workTest.get(0);
    Vector<T,U> work(lwork);
    dsyev_( "N", "L", &n, A.get(), &n, w.get(), work.get(), &lwork, &info );
    assert(info == 0);
}

template<typename T, typename U>
Vector<T,U> eigenValuesAndVectorsSymmetricMatrix(Matrix<T,U> &A) {
    Vector<T,U> w(A.getNx());
    U n = A.getNx();
    Vector<T,U> workTest(1);
    U lwork = -1, info = 0;
    
//     dsyev_( "V", "U", &n, A.get(), &n, w.get(), workTest.get(), &lwork, &info );
//     lwork= workTest.get(0);
//     Vector<T,U> work(lwork);
//     dsyev_( "V", "U", &n, A.get(), &n, w.get(), work.get(), &lwork, &info );
    
    Rsyev("V", "U", (mpackint)n, A.get(), (mpackint)n, w.get(), workTest.get(), (mpackint)lwork, &info);
    lwork = (int) workTest.get(0).get_d();
    Vector<T,U> work(lwork);
    Rsyev("V", "U", n, A.get(), n, w.get(), work.get(), lwork, &info);
    
    assert(info == 0);
    A.transposeInPlace();
    return w;
}

template<typename T, typename U>
std::vector<Matrix<T,U> > exponentiateMatrix(std::vector<T> &timeVec, T gamma, T h, U n, GridGenerator<T,U> &grid, T threshold) {
    BetaAndDeltaGridGenerator<T,U> constants(grid, gamma, h);
    
    std::auto_ptr<Matrix<T,U> > Q = computeSymmetricMatrixFromConstants(n, constants);
    Vector<T,U> w = eigenValuesAndVectorsSymmetricMatrix(*Q);
    Matrix<T,U> QT(*Q); QT.transposeInPlace();
    
    std::auto_ptr<Matrix<T,U> > DO = prod(constants.getSmallD(),*Q);
    std::auto_ptr<Matrix<T,U> > DOT = prod(QT, constants.getInvSmallD());
    
    std::auto_ptr<Matrix<T,U> > shouldBeId = prod(*DO,*DOT);
    std::vector<Matrix<T,U> > expVec;
    T diff = norm(*shouldBeId)-trace(*shouldBeId);
    if (!( diff <= threshold && diff >= -threshold)) {
        std::cout << "Norm(shouldBeId) - trace(shouldBeId) = " << diff << " , while threshold is = " << threshold << std::endl;
        return expVec;
    }
    
    std::auto_ptr<Matrix<T,U> > exponential = std::auto_ptr<Matrix<T,U> >(new Matrix<T,U>(w.getNx()));
    
    for (U iT = 0; iT < (U)timeVec.size(); ++iT) {
        T time = timeVec[iT];
        assert(time >= T());
        if (time == T()) {
            Matrix<T,U> id(Q->getNx()); id.setToConst(T());
            for (U iX = 0; iX < id.getNx(); ++iX) id.get(iX,iX) = (T)1;
            expVec.push_back(id);
        } else {
            exponential->setToConst(T());
            for (U iX = 0; iX < w.getNx(); ++iX) exponential->get(iX,iX) = exp(w.get(iX)*time);
            
            exponential = prod(*DO,*prod(*exponential,*DOT));
            if (!ensursePositivityAndOneSmallness(*exponential, threshold)) {
                std::cout << "Positiveness and smallness than one is not ensured, while threshold is = " << threshold << std::endl;
                expVec.clear();
                return expVec;
            }
            
            expVec.push_back(*exponential);
        }
    }
    
    return expVec;
}

template<typename T, typename U>
void computeExpTime(T time, Matrix<T,U> &exponential, Matrix<T,U> &bigV, Matrix<T,U> &lambda, Matrix<T,U> &lambdaPrime, Matrix<T,U> &id, Matrix<T,U> &DO, Matrix<T,U> &DOT) {
    exponential.setToConst(T());
    for (U iX = 0; iX < exponential.getNx(); ++iX) exponential.get(iX,iX) = exp(lambda.get(iX,iX)*time);
    
    std::auto_ptr<Matrix<T,U> > bigV_t   = prod(bigV,time);
    std::auto_ptr<Matrix<T,U> > lambda_t = prod(lambda,time);
    
    exponential = *prod(*prod(DO,*add(*add(*bigV_t, *prod(*prod(*sub(*sub(exponential, id),*lambda_t),lambdaPrime),bigV)),exponential)),DOT);
}

template<typename T, typename U>
std::vector<Matrix<T,U> > exponentiateBigMatrix(std::vector<T> &timeVec, T gamma, T h, U n, GridGenerator<T,U> &grid, T threshold) {
    BetaAndDeltaGridGenerator<T,U> constants(grid, gamma, h);
    
    std::auto_ptr<Matrix<T,U> > smallS = computeSymmetricMatrixFromConstants(n, constants);
    std::auto_ptr<Matrix<T,U> > bigR   = computeAlmostSymmetricMatrixFromConstants(n, constants);

    Vector<T,U> w                      = eigenValuesAndVectorsSymmetricMatrix(*smallS);
    std::auto_ptr<Matrix<T,U> > bigO   = computeBigO(*smallS);

    Matrix<T,U> bigOT(*bigO);
    bigOT.transposeInPlace();
    std::auto_ptr<Matrix<T,U> > bigT   = prod(bigOT, *prod(*bigR,*bigO));
    
    Matrix<T,U> lambda(n), lambdaPrime(n);
    lambda.setToConst(T()); lambdaPrime.setToConst(T());
    for (U iX = 0; iX < w.getNx(); ++iX) {
        lambda.get(iX+1,iX+1)      = w.get(iX);
        lambdaPrime.get(iX+1,iX+1) = (T)1/w.get(iX);
    }

    std::auto_ptr<Matrix<T,U> > bigV = sub(*bigT,lambda);
    
    Matrix<T,U> id(n); id.setToConst(T());
    for (U iX = 0; iX < id.getNx(); ++iX) id.get(iX,iX) = (T)1; // identity matrix
    
    std::auto_ptr<Matrix<T,U> > DO = prod(constants.getBigD(),*bigO);
    std::auto_ptr<Matrix<T,U> > DOT = prod(bigOT, constants.getInvBigD());
        
    std::auto_ptr<Matrix<T,U> > shouldBeId = prod(*DO,*DOT);
    std::vector<Matrix<T,U> > expVec;
    T diff = norm(*shouldBeId)-trace(*shouldBeId);
    if (!( diff <= threshold && diff >= -threshold)) {
        std::cout << "Norm(shouldBeId) - trace(shouldBeId) = " << diff << " , while threshold is = " << threshold << std::endl;
        return expVec;
    }

    Matrix<T,U> exponential(n);
//     T totTime = T();
    for (U iT = 0; iT < (U)timeVec.size(); ++iT) {
        T time = timeVec[iT];
//         totTime += time;
        assert(time >= T());
        if (time == T()) {
            Matrix<T,U> id(n); id.setToConst(T());
            for (U iX = 0; iX < id.getNx(); ++iX) id.get(iX,iX) = (T)1;
            expVec.push_back(id);
        } else {
            computeExpTime(time, exponential, *bigV, lambda, lambdaPrime, id, *DO, *DOT);
            if (!ensursePositivityAndOneSmallness(exponential, threshold)) {
                std::cout << "Positiveness and smallness than one is not ensured, while threshold is = " << threshold << std::endl;
                expVec.clear();
                return expVec;
            }
            
            expVec.push_back(exponential);
        }
    }
    
    return expVec;
}

#endif // EIGENVALUE_CALCULATIONS_HH
