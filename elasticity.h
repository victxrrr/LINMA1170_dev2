#include "matrix.h"

#ifndef _ELASTICITY_H_
#define _ELASTICITY_H_

void p1_stiffness_matrix_plane_stress (double E, double nu, double dphi[6], double det, double S[6][6]);
void p1_geometry(const double *x, double *detptr, double dxidx[2][2], double *dphi);
static void p1_stiffness_matrix (double dphi[6], double h[3][3] , double det, double S[6][6]);
static void hooke_plain_stress (double E, double nu, double h[3][3]);
void p1_stiffness_matrix_plane_stress (double E, double nu, double dphi[6], double det, double S[6][6]);
void p1_geometry(const double *x, double *detptr, double dxidx[2][2], double *dphi);
int get_triangles_p1 (size_t ** elementTags, size_t * elementTags_n, size_t ** nodeTags, size_t * nodeTags_n);
int get_nodes (size_t ** nodeTags, size_t * nodeTags_n,
		double ** coord, size_t * coord_n);
Matrix * compute_stiffness (
	int nNodes,
	double *coord,
	size_t ntriangles,
	size_t *triangleNodes);
void boundary_conditions (Matrix *K, double *RHS);
void assemble_system(Triplet** triplets, int* n_triplets, double** RHS, double** coord, size_t** gmsh_num, int* n_nodes);
void visualize_in_gmsh(double* SOL, size_t* gmsh_num, int n_nodes);

#endif
