#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <mblas_gmp.h>
#include <mlapack_gmp.h>

#include "headers.h"
#include "headers.hh"


using namespace std;

extern "C" void dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda,
                        double *w, double *work, int *lwork, int *info );

int main(int argc, char *argv[]) {
    
    if (argc != 5) {
        cout << "Error ! The input of the simulation should be : " << endl;
        cout << "1. gamma." << endl;
        cout << "2. h." << endl;
        cout << "3. H." << endl;
        cout << "4. t." << endl;
        exit(1);
    }
    
    double gamma = atof(argv[1]); 
    double h     = atof(argv[2]);
    int H        = atoi(argv[3]);
    
    double t = atof(argv[4]);
    Matrix<double,int> Q = exponentiateMatrix(t, gamma, h, H);
//     Q.print(6);
    for (int iX = 0; iX < Q.getNx(); ++iX) {
        for (int iY = 0; iY < Q.getNy(); ++iY) {
            if (/*Q.get(iX,iY) < double() || */Q.get(iX,iY) > (double)1) {
                cout << Q.get(iX,iY) << " " << endl;
            }
        }
    }
//     
//     cout << "Printing output Q." << endl;
//     Q.print();
//     cout << "Printing original R." << endl;
//     R.print();
// 
//     for (int iX = 0; iX < Q.getNx(); ++iX) {
//         Vector<double> vec(Q.col(iX));
//         for (int iY = 0; iY < Q.getNy(); ++iY) {
//             vec.get(iY) *= w.get(iX);
//         }
//         cout << "Printing rhs eigen vector " << iX << endl;
//         vec.print();
//     }
// 
//     for (int iX = 0; iX < Q.getNx(); ++iX) {
//         Vector<double> ini(R.row(iX)), res(n);
//         for (int iY = 0; iY < Q.getNy(); ++iY) {
//             Vector<double> ev(Q.col(iY));
//             res.get(iY) = 0.0;
//             for (int iD = 0; iD < n; ++iD) {
//                 res.get(iY) += ini.get(iD)*ev.get(iD);
//             }
//         }
//         cout << "Printing lhs eigen vector " << iX << endl;
//         res.print();
//     }
        
    
    return 0;
    
}
