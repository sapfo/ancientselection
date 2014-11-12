#ifndef CONSTANTS_HELPERS_H
#define CONSTANTS_HELPERS_H

#include "matrixOperations.h"
#include "matrixOperations.hh"

#include <cmath>

template<typename T, typename U>
class GridGenerator {
public :
    GridGenerator(U H_) : H(H_), xx(H)
    { }
    virtual ~GridGenerator() { };
    
public :
    const U getN() const {
        return H;
    }
    Vector<T,U> &get() {
        return xx;
    }
    
    Vector<T,U> const &get() const {
        return xx;
    }
    
    T &get(U iX) {
        assert(0 <= iX && iX < H);
        return xx.get(iX);
    }
    
    T const &get(U iX) const {
        assert(0 <= iX && iX < H);
        return xx.get(iX);
    }
private :
    U H;
    Vector<T,U> xx;
};

// Generation of the regular grid
template<typename T, typename U>
class RegularGridGenerator : public GridGenerator<T,U> {
public :
    RegularGridGenerator(U H_) : GridGenerator<T,U>(H_) {
        for (U iX = 0; iX < this->getN(); ++iX) {
            this->get(iX) = (T)iX / (T)(this->getN()-1);
        }
    }
};

// Generation of the quadratic grid
template<typename T, typename U>
class QuadraticGridGenerator : public GridGenerator<T,U> {
public :
    QuadraticGridGenerator(U H_, T Ne_) : GridGenerator<T,U>(H_), Ne(Ne_) {
        U smallPts = (U)(this->getN() / 10) + ((U(this->getN() / 10)+1) % 2);
        U largePts = this->getN() - smallPts+1;
        
        Vector<T,U> gridSmall(smallPts);
        for (U iX = 0; iX < smallPts; ++iX) {
            gridSmall.get(iX) = (T)iX / (T)(smallPts-1)*(T)1/Ne;
        }
//         U startFreq = (smallPts-1)/2;
        
        Vector<T,U> gridLarge(largePts);
        for (U iX = 0; iX < largePts; ++iX) {
            gridLarge.get(iX) = (T)iX / (T)(largePts-1);
        }
        T dq = gridLarge.get(1) - gridLarge.get(0);
        T xStart = gridSmall.get(gridSmall.getNx()-1);
        T dxStart = gridSmall.get(gridSmall.getNx()-1)-gridSmall.get(gridSmall.getNx()-2);
        
        for (U iN = 0; iN < smallPts-1; ++iN) {
            this->get(iN) = gridSmall.get(iN);
        }
        
        for (U iN = smallPts-1; iN < this->getN(); ++iN) {
            U iB = iN-smallPts+1;
            
            T q = gridLarge.get(iB);
            T b = -(T)3*(-dq + dxStart + dq*xStart)/dq;
            T a = -(T)2 * b/(T)3;
            
            this->get(iN) = a*q*q*q+b*q*q+dxStart/dq*q+xStart;
        }
    }
private :
    T Ne;
};

// beta and delta coefficients generator
template<typename T, typename U>
class BetaAndDeltaGridGenerator {
public :
    BetaAndDeltaGridGenerator(const GridGenerator<T,U> &xx, T gamma, T h );
public :
    // beta vector (r/w) cf. notations of AS
    Vector<T,U> &getBeta() {
        return beta;
    }
    
    Vector<T,U> &getDelta() {
        return delta;
    }
    
    Vector<T,U> &getEta() {
        return eta;
    }
    
    Matrix<T,U> &getSmallD() {
        return smallD;
    }
    
    Matrix<T,U> &getInvSmallD() {
        return invSmallD;
    }
    
    Matrix<T,U> &getBigD() {
        return bigD;
    }
    
    Matrix<T,U> &getInvBigD() {
        return invBigD;
    }
    
private :
    U H;
    Vector<T,U> beta, delta, eta;
    Matrix<T,U> smallD, invSmallD, bigD, invBigD;
};

template<typename T, typename U>
BetaAndDeltaGridGenerator<T,U>::BetaAndDeltaGridGenerator(const GridGenerator<T,U> &xx, T gamma, T h ) 
        : H(xx.getN()), beta(H-2), delta(H-2), eta(H-2), smallD(H-2), invSmallD(H-2), bigD(H), invBigD(H) {
    
    for (U iX = 1; iX < H-1; ++iX) {
        T xup   = xx.get(iX+1);
        T x     = xx.get(iX);
        T xdown = xx.get(iX-1);
        
        beta.get(iX-1)  =  ((-(T)1 + x) * x * (-(T)1-(x*x)*gamma + h*(-(T)1 + (T)2*x)*(x - xdown)*gamma+x*xdown*gamma))
                            / ((x - xup)*(xdown -xup));
        delta.get(iX-1) = -((-(T)1 + x)*x*(-(T)1-(x*x)*gamma + h*(-(T)1 + (T)2*x)*(x - xup)*gamma+x*xup*gamma ))
                            / ((x - xdown)*(xdown - xup));
                            
        eta.get(iX-1) = -(beta.get(iX-1) + delta.get(iX-1));
    }
    
    smallD.setToConst(T()); invSmallD.setToConst(T());
    smallD.get(0,0) = (T)1;
    invSmallD.get(0,0) = smallD.get(0,0);
    for (U iX = 1; iX < H-2; ++iX) {
        T div = delta.get(iX)/beta.get(iX-1);
        assert(delta.get(iX) >= T() && beta.get(iX-1) >= T());
        smallD.get(iX,iX) = sqrt(div)*smallD.get(iX-1,iX-1);
        invSmallD.get(iX,iX) = (T)1/smallD.get(iX,iX);
    }
    
    bigD.setToConst(T()); invBigD.setToConst(T());
    bigD.get(0,0) = (T)1; bigD.get(H-1,H-1) = (T)1;
    invBigD.get(0,0) = (T)1; invBigD.get(H-1,H-1) = (T)1;
    for (U iX = 1; iX < H-1; ++iX) {
        bigD.get(iX,iX) = smallD.get(iX-1,iX-1);
        invBigD.get(iX,iX) = (T)1 / bigD.get(iX,iX);
    }
}

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> >  computeSymmetricMatrixFromConstants(U size, BetaAndDeltaGridGenerator<T,U> &constants) {
    Vector<T,U> beta = constants.getBeta();
    Vector<T,U> delta = constants.getDelta();
    Vector<T,U> eta = constants.getEta();
    
    Matrix<T,U> q(size-2);
    q.setToConst(T());
    // upper part of the matrix filling
    for (U iX = 0; iX < size-3; ++iX) {
        q.get(iX,iX+1) = beta.get(iX);
    }
    // lower part of the matrix filling
    for (U iX = 0; iX < size-3; ++iX) {
        q.get(iX+1,iX) = delta.get(iX+1);
    }
    // diagonal of the matrix filling
    for (U iX = 0; iX < size-2; ++iX) {
        q.get(iX,iX) = eta.get(iX);
    }
    
    return prod(constants.getInvSmallD(),*prod(q,constants.getSmallD()));
}

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > computeAlmostSymmetricMatrixFromConstants(U size, BetaAndDeltaGridGenerator<T,U> &constants) {
    Vector<T,U> beta = constants.getBeta();
    Vector<T,U> delta = constants.getDelta();
    Vector<T,U> eta = constants.getEta();
    
    Matrix<T,U> bigQ(size);
    bigQ.setToConst(T());
    for (U iX = 1; iX < size-1; ++iX) {
        bigQ.get(iX,iX+1) = beta.get(iX-1);
        bigQ.get(iX,iX-1) = delta.get(iX-1);
        bigQ.get(iX,iX)   = eta.get(iX-1);
    }
    bigQ.print("bigQ.dat");
    
    return prod(constants.getInvBigD(), *prod(bigQ, constants.getBigD()));
}

template<typename T, typename U>
std::auto_ptr<Matrix<T,U> > computeBigO(const Matrix<T,U> &smallO) {
    Matrix<T,U> *bigO = new Matrix<T,U>(smallO.getNx()+2, smallO.getNy()+2);
    bigO->setToConst(T());
    bigO->get(0,0) = (T)1;
    bigO->get(bigO->getNx()-1,bigO->getNy()-1) = (T)1;
    for (U iX = 0; iX < smallO.getNx(); ++iX) {
        for (U iY = 0; iY < smallO.getNy(); ++iY) {
            bigO->get(iX+1,iY+1) = smallO.get(iX,iY);
        }
    }
    return std::auto_ptr<Matrix<T,U> >(bigO);
}

#endif // CONSTANTS_HELPERS_H
