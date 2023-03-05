#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "matrix.h"


/**
 * Alloue une matrice de dimension m x n
 * @param m: nombre de lignes
 * @param n: nombre de colonnes
 * @return pointeur vers une nouvelle structure de données allouée Matrix
*/
Matrix * allocate_matrix(int m, int n) {
	Matrix *mat = (Matrix *) malloc(sizeof(Matrix));
	if (mat == NULL) {
		strerror(errno);
	}
	mat->m = m; mat->n = n;

	double *data = (double *) malloc(sizeof(double)*m*n);
	double **a = (double **) malloc(sizeof(double *)*m);
	if (data == NULL || a == NULL) {
		strerror(errno);
	}
	for (int i = 0; i < m; i++) *(a + i) = data + i*n;
	mat->data = data;
	mat->a = a;
	return mat;
}

/**
 * Alloue une matrice bande de dimension m x m et largeur de bande k
 * @param m: nombre de lignes et de colonnes
 * @param k: largeur de bande
 * @return pointeur vers une nouvelle structure de données allouée BandMatrix
*/
BandMatrix * allocate_band_matrix(int m, int k) {
    BandMatrix *mat = (BandMatrix *) malloc(sizeof(BandMatrix));
    if (mat == NULL) {
        strerror(errno);
    }
    mat->m = m; mat->k = k;

    double *data = (double *) malloc(sizeof(double)*m*(2*k+1));
    double **a = (double **) malloc(sizeof(double *)*m);
    if (data == NULL || a == NULL) {
        strerror(errno);
    }
    for (int i = 0; i < m; i++) *(a + i) = data + i*(2*k+1);
    mat->data = data;
    mat->a = a;
    return mat;
}

/**
 * Libère du heap une structure Matrix précédemment allouée
 * @param mat: structure à libérer
*/
void free_matrix(Matrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}

/**
 * Libère du heap une structure BandMatrix précédemment allouée
 * @param mat: structure à libérer
*/
void free_band_matrix(BandMatrix * mat) {
    free(mat->a);
    free(mat->data);
    free(mat);
}

/**
 * Calcule une permutation pour réduire le fill-in.
 * @param perm: tableau de taille n_nodes contenant la permutation à calculer
 * @param coord: tableau de taille 2*n_nodes contenant les coordonnées des noeuds
 * @param n_nodes: nombre de noeuds
 * @param triplets: tableau de taille n_triplets contenant les triplets des NNZ de la matrice
 * @param n_triplets: nombre de triplets
 * @return nombre de permutations effectuées
*/
int compute_permutation(int * perm, double * coord, int n_nodes, Triplet * triplets, int n_triplets) {
    return 0;
}


/**
 * Imprime le contenu du vecteur v sur la sortie standard
 * @param v: vecteur de taille n
 * @param n: taille du vecteur
*/
void print_vector(double * v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%.3f\t", v[i]);
    }
    printf("\n");
}

/**
 * Imprime le contenu de la matrice A sur la sortie standard
 * @param A: matrice de dimension m x n
*/
void print_matrix(Matrix * A) {

	for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->n; j++) {
            printf("%.3f\t", A->a[i][j]);
        }
        printf("\n");
    }
	printf("\n");
}