#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <mblas_gmp.h>
#include <mlapack_gmp.h>
#include <memory>
#include <omp.h>

#include "headers.h"
#include "headers.hh"

using namespace std;

typedef mpf_class T;
typedef mpackint U;

template<typename Ty>
std::string stringify(Ty number) {
    std::ostringstream oss;
    if (!(oss << number)) {
        std::cout << "Error in function stringify" << std::endl;
        exit(1);
    }

    return oss.str();
}

int main(int argc, char *argv[]) {
    
    if (argc < 7) {
        cout << "Error ! The input of the simulation should be : " << endl;
        cout << "1. number precision." << endl;
        cout << "2. threshold." << endl;
        cout << "3. gamma." << endl;
        cout << "4. h." << endl;
        cout << "5. H." << endl;
        cout << "6. if quadratic grid is use (0 no, 1 yes)." << endl;
        cout << "7. if 6 == 0 begining of times, else Ne value." << endl;
        cout << "8. If 6 == 1 begining of times." << endl;
        exit(1);
    }
    
    int prec          = atoi(argv[1]);
    //initialization of GMP
    mpf_set_default_prec(prec);
    
    mpf_class epsilon = atof(argv[2]);
    mpf_class gamma   = atof(argv[3]); 
    mpf_class h       = atof(argv[4]);
    mpackint H        = atoi(argv[5]);
    int quadratic     = atoi(argv[6]);
    
    GridGenerator<T,U> *grid;
    if (quadratic == 0) {
        grid = new RegularGridGenerator<T,U>(H);
    } else if (quadratic == 1) {
        mpackint Ne     = atoi(argv[7]);
        grid = new QuadraticGridGenerator<T,U>(H, Ne);
    } else {
        std::cout << "Invalid fourth argument. Must be 0 or one." << std::endl;
        exit(1);
    }
    
    std::vector<T> time;
    for (U iT = 7+quadratic; iT < argc; ++iT) {
        T timeTmp = atof(argv[iT]);
        time.push_back(timeTmp);
    }
    
    std::vector<Matrix<T,U> > expBigQ = exponentiateMatrix(time, gamma, h, H, *grid, epsilon);

    for (U iT = 0; iT < (U)expBigQ.size(); ++iT) {
        expBigQ[iT].print("exp_Qt_"+stringify<U>(iT)+".dat");
    }
    
    return 0;
    
}
