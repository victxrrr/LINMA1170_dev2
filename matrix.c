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
 * Alloue une matrice symétrique bande de dimension m x m et largeur de bande k
 * @param m: nombre de lignes et de colonnes
 * @param k: largeur de bande
 * @return pointeur vers une nouvelle structure de données allouée SymBandMatrix
*/
SymBandMatrix * allocate_sym_band_matrix(int m, int k) {
    SymBandMatrix *mat = (SymBandMatrix *) malloc(sizeof(SymBandMatrix));
    if (mat == NULL) {
        strerror(errno);
    }
    mat->m = m; mat->k = k;

    double *data = (double *) malloc(sizeof(double)*m*(k+1));
    double **a = (double **) malloc(sizeof(double *)*m);
    if (data == NULL || a == NULL) {
        strerror(errno);
    }
    for (int i = 0; i < m; i++) *(a + i) = data + k*i;
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
 * Libère du heap une structure SymBandMatrix précédemment allouée
 * @param mat: structure à libérer
*/
void free_sym_band_matrix(SymBandMatrix * mat) {
    free(mat->a);
    free(mat->data);
    free(mat);
}

/**
 * Renvoie le minimum de deux entiers
*/
int min(int a, int b) {
    return a < b ? a : b;
}

/**
 * Renvoie le maximum de deux entiers
*/
int max(int a, int b) {
    return a > b ? a : b;
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
    for (int i = 0; i < n_nodes; i++) {
        perm[2*nodes[i].index] = 2*i;
        perm[2*nodes[i].index+1] = 2*i+1;
    }

    // update triplets
    for (int i = 0; i < n_triplets; i++) {
        triplets[i].i = perm[triplets[i].i];
        triplets[i].j = perm[triplets[i].j];
    }

    free(nodes);
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
 * Imprime le contenu de la matrice bande A sur la sortie standard
 * @param A: matrice bande de dimension m x m et largeur de bande k
*/
void print_band_matrix(BandMatrix * A) {
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (abs(i - j) <= A->k) printf("%.2e\t", A->a[i][j]);
            else printf("%.2e\t", 0.0);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * Imprime le contenu de la matrice bande symétrique A sur la sortie standard
 * @param A: matrice bande symétrique de dimension m x m et largeur de bande k
*/
void print_sym_band_matrix(SymBandMatrix * A) {
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (i <= j && abs(i - j) <= A->k) printf("%.2e\t", A->a[i][j]);
            else printf("%.2e\t", 0.0);
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

/**
 * Ecrit le contenu d'une matrice dans un fichier CSV
 * @param A: matrice de dimension m x n
 * @param filename: nom du fichier
*/
void matrix_to_csv(Matrix * A, char * filename) {
    FILE *fp;
    fp = fopen(filename, "w");
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->n; j++) {
            if (j==A->n-1) fprintf(fp, "%lf", A->a[i][j]);
     		else fprintf(fp, "%lf, ", A->a[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/**
 * Ecrit le contenu d'une matrice bande dans un fichier CSV
 * @param A: matrice carrée de dimension m et de largeur de bande k
 * @param filename: nom du fichier
*/
void bandmatrix_to_csv(BandMatrix * A, char * filename) {
    FILE *fp;
    fp = fopen(filename, "w");
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (abs(j - i) > A->k) {
                if (j == A->m - 1) fprintf(fp, "%lf", 0.0);
                else fprintf(fp, "%lf,", 0.0);
            } else {
                if (j == A->m - 1) fprintf(fp, "%lf", A->a[i][j]);
                else fprintf(fp, "%lf,", A->a[i][j]);
            }
        }
        fprintf(fp, "\n");
    }
}

void symbandmatrix_to_csv(SymBandMatrix * A, char * filename) {
    FILE *fp;
    fp = fopen(filename, "w");
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (abs(i - j) <= A->k) {
                double aij = (i <= j) ? A->a[i][j] : 0.0;
                if (j == A->m - 1) fprintf(fp, "%lf", aij);
                else fprintf(fp, "%lf,", aij);
            } else {
                if (j == A->m - 1) fprintf(fp, "%lf", 0.0);
                else fprintf(fp, "%lf,", 0.0);
            }
        }
        fprintf(fp, "\n");
    }
}

/**
 * Fonction auxiliaire calculant le produit matriciel C = A@B
 * @param A: matrice de dimension m x p
 * @param B: matrice de dimension p x n
 * @param C: matrice de dimension m x n
 * @return -1 si les dimensions sont incompatibles, 0 sinon
*/
int mult_matrix(Matrix * A, Matrix * B, Matrix * C) {
	int m = A->m, p = A->n, n = B->n;
    if (B->m != p) return -1;
    if (C->m != m || C->n != n) return -1;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C->data[i * n + j] = 0.0;
			for (int k = 0; k < p; k++) {
				C->data[i * n + j] += A->a[i][k] * B->a[k][j];
			}
		}
	}
	return 0;
}