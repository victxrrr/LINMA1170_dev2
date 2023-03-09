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
    for (int i = 0; i < m; i++) *(a + i) = data + k*(2*i+1);
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
 * Imprime le contenu d'une matrice sur la sortie standard
 * @param A: matrice à imprimer
*/
int cmpfunc (const void * a, const void * b) {
    Node *aNode = (Node *) a;
    Node *bNode = (Node *) b;
    if (aNode->x > bNode->x) {
        return 1;
    } else if (aNode->x < bNode->x) {
        return -1;
    } else {
        if (aNode->y > bNode->y) {
            return 1;
        } else if (aNode->y < bNode->y) {
            return -1;
        } else {
            return 0;
        }
    }
}

/**
 * Calcule une permutation pour réduire le fill-in.
 * La stratégie de permutation consiste à trier les noeuds par ordre croissant de coordonnées.
 * @param perm: tableau de taille 2*n_nodes contenant la permutation à calculer
 * @param coord: tableau de taille 2*n_nodes contenant les coordonnées des noeuds
 * @param n_nodes: nombre de noeuds
 * @param triplets: tableau de taille n_triplets contenant les triplets des NNZ de la matrice
 * @param n_triplets: nombre de triplets
 * @return nombre de permutations effectuées
*/
int compute_permutation(int * perm, double * coord, int n_nodes, Triplet * triplets, int n_triplets) {

    Node * nodes = (Node *) malloc(sizeof(Node)*n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        nodes[i].x = coord[2*i];
        nodes[i].y = coord[2*i+1];
        nodes[i].index = i;
    }
    qsort(nodes, n_nodes, sizeof(Node), cmpfunc);

    // save the permutation
    int *perm_helper = (int *) malloc(sizeof(int)*2*n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        perm[2*nodes[i].index] = 2*i;
        perm[2*nodes[i].index+1] = 2*i+1;
    }

    // update triplets
    for (int i = 0; i < n_triplets; i++) {
        triplets[i].i = perm[triplets[i].i];
        triplets[i].j = perm[triplets[i].j];
    }

    return 0;
}


/**
 * Imprime le contenu du vecteur v sur la sortie standard
 * @param v: vecteur de taille n
 * @param n: taille du vecteur
*/
void print_vector(double * v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%.2e\t", v[i]);
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
            printf("%.2e\t", A->a[i][j]);
        }
        printf("\n");
    }
	printf("\n");
}

/**
 * Imprime le contenu d'un tableau de triplets sur la sortie standard
 * @param t: tableau de taille n
 * @param n: taille du tableau
*/
void print_triplets(Triplet * t, int n) {
    for (int i = 0; i < n; i++) {
        printf("%d\t%d\t%.2e\n", t[i].i, t[i].j, t[i].val);
    }
    printf("\n");
}

// int main(){

//     // BandMatrix *A = allocate_band_matrix(4, 1);

//     // double tab[4][4] = {{1.0, 2.0, 0.0, 0.0},
//     //                     {2.0, 3.0, 4.0, 0.0},
//     //                     {0.0, 4.0, 5.0, 6.0},
//     //                     {0.0, 0.0, 6.0, 7.0}};
//     // for (int i = 0; i < 4; i++) {
//     //     for (int j = 0; j < 4; j++) {
//     //         if (tab[i][j] != 0.0) A->a[i][j] = tab[i][j];
//     //     }
//     // }
//     // print_vector(A->data, 12);

//     // for (int i = 0; i < 4; i++) {
//     //     for (int j = 0; j < 4; j++) {
//     //         printf("%f\t", A->a[i][j]);
//     //     }
//     //     printf("\n");
//     // }

//     // free_band_matrix(A);

//     int N = 5, k = 2;
//     BandMatrix *A = allocate_band_matrix(N, k);
//     double tab[5][5] = {{1.0, 2.0, 4.0, 0.0, 0.0},
//                         {2.0, 3.0, 4.0, 0.0, 0.0},
//                         {4.0, 4.0, 5.0, 6.0, 1.0},
//                         {0.0, 3.0, 6.0, 7.0, 2.0},
//                         {0.0, 0.0, 1.0, 2.0, 3.0}};
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             if (tab[i][j] != 0.0) A->a[i][j] = tab[i][j];
//         }
//     }
//     print_vector(A->data, 12);

//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             printf("%f\t", A->a[i][j]);
//         }
//         printf("\n");
//     }

//     free_band_matrix(A);


//     return 0;
// }