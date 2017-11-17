/*******************************************************
#
# Panwen Wang, April 29th, 2009
# pwwang AT pwwang.com
# get eigenvalue and eigenvector of
# a real symmetric matrix like
#  1  -1  0
# -1   1  2
#  0   2  1
#
# Algorthm:
# 1. get a nonzero and non-diagonal element a[i,j] of matrix A, usually
#    the maximun absolute value of non-diagonal elements of A (max(A));
# 2. give sin(phi) and cos(phi) by
#    ( a[j,j]-a[i,i] )*sin(2*phi) + 2*a[i,j]*cos(2*phi) = 0 ;
#    => tan(2*phi) = 2*a[i,j]/(a[i,i]-a[j,j])
#    => phi = arctan(2*a[i,j]/(a[i,i]-a[j,j]))/2
# 3. get elements of a new matrix A1(a1[i,j]) by
#    a1[i,i] = a[i,i]*cos^2(phi) + a[j,j]*sin^2(phi) + 2*a[i,j]*cos(phi)*sin(phi)
#    a1[j,j] = a[i,i]*sin^2(phi) + a[j,j]*cos^2(phi) - 2*a[i,j]*cos(phi)*sin(phi)
#    a1[i,l] = a1[l,i] = a[i,l]*cos(phi) + a[j,l]*sin(phi)       ( l!=i,j )
#    a1[j,l] = a1[l,j] = -a[i,l]*sin(phi) + a[j,l]*cos(phi)      ( l!=i,j )
#    a1[l,m] = a1[m,l] = a[m,l]                                  ( m,l != i,j )
#    a1[i,j] = a1[j,i] = { (a[j,j]-a[i,i])*sin(2*phi) }/2 + a[i,j]*(cos^2(phi) - sin^2(phi))
# 4. let A1 be the substitution of A, repeat step 1,2,3 and get A2, and A3,A4,...,An can be
#    obtained by the same way. Calculation ceases if max(An) is less than the given threshold.
#
# Usage:
#        JACOBI_ARRAY ja(3), value;
#        JACOBI_MATRIX jm, vector;
#        ja[0] = 1; ja[1] = -1; ja[2] = 0;
#        jm.push_back(ja);
#        ja[0] = -1; ja[1] = 1; ja[2] = 2;
#        jm.push_back(ja);
#        ja[0] = 0; ja[1] = 2; ja[0] = 1;
#        jm.push_back(ja);
#        Jacobi j;
#        j.setMatrix(jm);
#         value = j.getEigenvalues();
#        vector = j.getEigenvectors();
#
********************************************************/
#ifndef _JACOBI_H_
#define _JACOBI_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#define JACOBI_ARRAY vector<float>
#define JACOBI_MATRIX vector< vector<float> >
#define JACOBI_vint vector<float>::size_type

using namespace std;

class Jacobi{
public:
    JACOBI_MATRIX matrix;
    JACOBI_ARRAY  eigenvalues;
    JACOBI_MATRIX eigenvectors;
    float e;   //precision
    int dimen;
    bool debug;

    float max(const JACOBI_MATRIX &, int &, int &);
    JACOBI_MATRIX newMatrix(const JACOBI_MATRIX &, const int &, const int &, const float &);
    void jacobi();

public:
    Jacobi();
    void setMatrix(const JACOBI_MATRIX &);
    void setE(const float &);
    void setDebug(const bool &d = true);
    JACOBI_ARRAY getEigenvalues();
    JACOBI_MATRIX getEigenvectors();
    void printEigen();
    void printMatrix();
    void printMatrix(const JACOBI_MATRIX &);
    void demo();
};

Jacobi::Jacobi():e(1e-5),debug(false){}

void Jacobi::setE(const float &ee){ e = ee; }

void Jacobi::setDebug(const bool &d){ debug = d; }

void Jacobi::setMatrix(const JACOBI_MATRIX & jm){
    dimen = jm.size();
    matrix = jm;
    eigenvalues.resize(dimen);
    eigenvectors.resize(dimen);
    //matrix.resize(dimen);
    //for(JACOBI_vint i=0; i<dimen; i++){ matrix[i] = jm[i];}
    jacobi();
}

JACOBI_ARRAY Jacobi::getEigenvalues(){return eigenvalues;}

JACOBI_MATRIX Jacobi::getEigenvectors(){return eigenvectors;}

float Jacobi::max(const JACOBI_MATRIX &jm, int &maxi, int &maxj){
    float m = jm[0][1];
    maxi=0, maxj=1;
    for(JACOBI_vint i=0; i<dimen; i++)
        for(JACOBI_vint j=0; j<dimen; j++)
            if(i!=j)
                if(fabs(jm[i][j])>m){
                    m = fabs(jm[i][j]);
                    maxi = i; maxj = j;
                }
    return m;
}

JACOBI_MATRIX Jacobi::newMatrix(const JACOBI_MATRIX &jm, const int &i, const int &j, const float &phi){
    JACOBI_MATRIX jm1(dimen);
    float sp = sin(phi), cp = cos(phi);
    for(int ii=0; ii<dimen; ii++){
        jm1[ii].resize(dimen);
        for(int jj=0; jj<dimen; jj++){
            if( ii==i ){     // row i
                if( jj==i ) jm1[ii][jj] = jm[i][i]*cp*cp + jm[j][j]*sp*sp + 2*jm[i][j]*cp*sp;
                else if( jj==j ) jm1[ii][jj] = (jm[j][j]-jm[i][i])*sp*cp + jm[i][j]*(cp*cp-sp*sp);
                else jm1[ii][jj] = jm[i][jj]*cp + jm[j][jj]*sp;
            } else if ( ii==j ) {// row j
                if( jj==i ) jm1[ii][jj] = (jm[j][j]-jm[i][i])*sp*cp + jm[i][j]*(cp*cp-sp*sp);
                else if( jj==j ) jm1[ii][jj] = jm[i][i]*sp*sp + jm[j][j]*cp*cp - 2*jm[i][j]*cp*sp;
                else jm1[ii][jj] = jm[j][jj]*cp - jm[i][jj]*sp;
            } else {            // row l ( l!=i,j )
                if( jj==i ) jm1[ii][jj] = jm[i][ii]*cp + jm[j][ii]*sp;
                else if( jj==j ) jm1[ii][jj] = jm[j][ii]*cp - jm[i][ii]*sp;
                else jm1[ii][jj] = jm[ii][jj];
            }
        }
    }
    return jm1;
}

void Jacobi::jacobi(){
    JACOBI_MATRIX jm = matrix;  // initial matrix
    float phi, m;
    int i,j;
    JACOBI_MATRIX tempVectors(dimen), vectors(dimen), cmat(dimen);
    for(JACOBI_vint x=0; x<dimen; x++){
        tempVectors[x].resize(dimen);
        for(JACOBI_vint y=0; y<dimen; y++){
            if(x==y) { tempVectors[x][y] = 1.0; }
            else { tempVectors[x][y] = 0.0; }
        }
    }

    for( ;; ){ // step 4
        m =  max(jm,i,j);  // step 1
        if(debug) printf("%-3d, %-3d, %.9f\n",i,j,m);
        if( m < e ) break;
        // if the maximum value of the matrix LE(Less than or Equal to) the limit, break
        phi = atan2(2*jm[i][j], jm[i][i] - jm[j][j]) / 2;  // step 2
        jm = newMatrix(jm,i,j,phi); // step 3
        for(JACOBI_vint x=0; x<dimen; x++){
            cmat[x].resize(dimen);
            vectors[x].resize(dimen);
            for(JACOBI_vint y=0; y<dimen; y++){
                if(x==y) { cmat[x][y] = 1.0; }
                else { cmat[x][y] = 0.0; }
                vectors[x][y] = 0.0;
            }
        }
        cmat[i][i] = cos(phi);
        cmat[j][j] = cos(phi);
        cmat[i][j] = -sin(phi);
        cmat[j][i] = sin(phi);
        for(JACOBI_vint x=0; x<dimen; x++){ // for eigenvectors
            for(JACOBI_vint y=0; y<dimen; y++){
                for(JACOBI_vint z=0; z<dimen; z++){
                    vectors[x][y] = vectors[x][y] + tempVectors[x][z] * cmat[z][y];
                }
                //print vectors[x][y],"-",tempVectors[x][y]," ";
            }
        }
        tempVectors = vectors;
        //@tempVectors = @vectors;
    }
    // the digonal elements are the eigenvalues
    for(JACOBI_vint x=0; x<dimen; x++) eigenvalues[x] = jm[x][x] ;

    for(JACOBI_vint x=0; x<dimen; x++){
        eigenvectors[x].resize(dimen);
        for(JACOBI_vint y=0; y<dimen; y++){
            eigenvectors[x][y] = vectors[y][x];
        }
    }
    if(debug) printEigen() ;
}

void Jacobi::printEigen(){
    for(JACOBI_vint i=0; i<dimen; i++){
        cout << "Eigenvalue [" << i << "]: " << eigenvalues[i] << "\n";
        cout << "Eigenvector[" << i << "]: [ ";
        for(JACOBI_vint j=0; j<dimen; j++){
            printf("%9.6f, ",eigenvectors[i][j]);
        }
        cout << " ] \n";
    }
}

void Jacobi::printMatrix(){ printMatrix(matrix); }

void Jacobi::printMatrix(const JACOBI_MATRIX &jm){
    for(JACOBI_vint i=0; i<dimen; i++){
        for(JACOBI_vint j=0; j<dimen; j++){
            printf("%9.6f, ", jm[i][j]);
        }
        cout << "\n";
    }
}

void Jacobi::demo(){
    dimen = 3;
    eigenvalues.resize(dimen);
    eigenvectors.resize(dimen);
    matrix.resize(dimen);
    debug = false;
    e = 1e-5;
    JACOBI_ARRAY ja(3);
    ja[0] = 2; ja[1] = -1; ja[2] = 0;
    matrix[0] = ja;
    ja[0] = -1; ja[1] = 2; ja[2] = -1;
    matrix[1] = ja;
    ja[0] = 0; ja[1] = -1; ja[2] = 2;
    matrix[2] = ja;
    printMatrix();
    jacobi();
    printEigen();
    /*
     2.000000, -1.000000,  0.000000,
    -1.000000,  2.000000, -1.000000,
     0.000000, -1.000000,  2.000000,
    Eigenvalue [0]: 3.41421
    Eigenvector[0]: [  0.500000, -0.707107,  0.500000,  ]
    Eigenvalue [1]: 2
    Eigenvector[1]: [  0.707107,  0.000000, -0.707107,  ]
    Eigenvalue [2]: 0.585786
    Eigenvector[2]: [  0.500000,  0.707107,  0.500000,  ]
    */
}

#endif
