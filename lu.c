#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lu.h"


/**
 * Décompose la matrice A en sa factorisation LU : A = L@U
 * @param A: matrice n x n à décomposer
 * @return -1 si l'algorithme rencontre un pivot nul, 0 sinon
*/
int lu(Matrix * A) {

    int n = A->m;

    for (int k = 0; k < n; k++) {
        double akk = A->a[k][k];
        if (fabs(akk) < EPS) return -1;
        for (int i = k + 1; i < n; i++) {
            A->a[i][k] /= akk;
            for (int j = k + 1; j < n; j++) {
                A->a[i][j] -= A->a[i][k]*A->a[k][j];
            }
        }
    }
    return 0;
}

/**
 * Résout LUx = y
 * @param LU: factorisation LU d'une certaine matrice
 * @param y: vecteur modifié contenant la solution x
 * @return 0 si la fonction réussit
*/
int solve(Matrix * LU, double * y){

    int n = LU->m;

    // forward substitution 
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < k; i++) {
            y[k] -= LU->a[k][i]*y[i];
        }
    }

    // backward substitution
    for (int k = n - 1; k >= 0; k--) {
        for (int i = k + 1; i < n; i++) {
            y[k] -= LU->a[k][i]*y[i];
        }
        y[k] /= LU->a[k][k];
    }

    return 0;
}

/**
 * Décompose la matrice A en sa factorisation de Choleski : A = R*@R
 * R est stocké dans A
 * @param A: matrice n x n définie positive et symétrique
 * @return 0
*/
int cholesky(Matrix * A) {

    int n = A->m;

    for (int k = 0; k < n; k++) {
        double Rkk = A->a[k][k];
        if (Rkk <= 0) return -1; // si A n'est pas définie positive
        for (int i = k + 1; i < n; i++) {
            double K = A->a[k][i]/Rkk;
            for (int j = i; j < n; j++) {
                A->a[i][j] -= A->a[k][j]*K;
            }
        }
        for (int i = k; i < n; i++) {
            A->a[k][i] = A->a[k][i]/sqrt(Rkk);
        }
    }
    return 0; 
}

/*  
 * Décompose la matrice bande A en sa factorisation LU : A = L@U
 * @param A: matrice m x m à décomposer
*/
int lu_band(BandMatrix * A) {

    int n = A->m, kmax = A->k;

    for (int k = 0; k < n; k++) {
        double akk = A->a[k][k];
        if (fabs(akk) < EPS) return -1;
        for (int i = k + 1; i <= min(k + kmax, n - 1); i++) {
            A->a[i][k] /= akk;
            for (int j = k + 1; j <= min(k + kmax, n - 1); j++) {
                A->a[i][j] -= A->a[i][k]*A->a[k][j];
            }
        }
    }
    return 0;
}

/*
* Résout LUx = y pour une matrice bande LU
* @param LU: factorisation LU d'une certaine matrice bande
* @param y: vecteur modifié contenant la solution x
*/
int solve_band(BandMatrix * LU, double * y){

    int n = LU->m, kmax = LU->k;

    // forward substitution 
    for (int k = 0; k < n; k++) {
        for (int i = max(0, k - kmax); i < k; i++) {
            y[k] -= LU->a[k][i]*y[i];
        }
    }

    // backward substitution
    for (int k = n - 1; k >= 0; k--) {
        for (int i = k + 1; i <= min(k + kmax, n - 1); i++) {
            y[k] -= LU->a[k][i]*y[i];
        }
        y[k] /= LU->a[k][k];
    }

    return 0;
}