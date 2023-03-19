#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../gmsh-sdk/include/gmshc.h"
#include "matrix.h"
#include "elasticity.h"
#include "lu.h"

int main (int argc, char *argv[]) {

	if (argc != 3) {
		printf("Usage: \n"
			"./solve_deformation <geo_file.geo> <meshSizeFactor>\n" 
			"---------------------------- \n\n"
			"- Use this .geo file : \n"
			"      square.geo \n"
			"- meshSizeFactor sets a size factor on the mesh size; for the square: \n "
			"      0. < meshSizeFactor < 1.\n \n");
		return -1;
	}

	int ierr; // Gmsh error code

	// Initialize Gmsh and load geometry
	gmshInitialize(argc, argv, 0, 0, &ierr);
	gmshOpen(argv[1], &ierr);

	// Set mesh size factor and generate mesh
	double meshSizeFactor;
	sscanf(argv[2],"%lf",&meshSizeFactor);
	gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
	gmshModelMeshGenerate(2, &ierr);

	// Setup linear system
	int n_nodes, n_triplets;
	size_t *gmsh_num;  // For final renumbering in Gmsh
	double *coord;     // Vector :  Coordinates of nodes (size 2*n_nodes)
	double *RHS;       // Vector : right-hand side of equation (size 2*n_nodes)
	Triplet* triplets; // Triplets of values to insert in matrix; each entry is a struct containing (i, j, val)
	assemble_system(&triplets, &n_triplets, &RHS, &coord, &gmsh_num, &n_nodes);

	// Build matrix
	// --- 1 ---
	Matrix * K = allocate_matrix(2*n_nodes, 2*n_nodes);
	for (int i = 0; i < K->m*K->n; i++) K->data[i] = 0.0;
	for (int t=0; t<n_triplets; t++) K->a[triplets[t].i][triplets[t].j] += triplets[t].val;

	// --- 2 ---
	int perm[2*n_nodes];
	compute_permutation(perm, coord, n_nodes, triplets, n_triplets);

	Matrix *permK = allocate_matrix(2*n_nodes, 2*n_nodes);
	for (int i = 0; i < permK->m*permK->n; i++) permK->data[i] = 0.0;

	int k = 0; // compute bandwidth
	for (int t = 0; t < n_triplets; t++) {
		if (abs(triplets[t].j - triplets[t].i) > k) k = abs(triplets[t].j - triplets[t].i);
		permK->a[triplets[t].i][triplets[t].j] += triplets[t].val;
	}

	double permRHS[2*n_nodes];
	for (int i = 0; i < 2*n_nodes; i++) permRHS[perm[i]] = RHS[i]; // permute RHS

	// --- 3 ---
	BandMatrix *bandK = allocate_band_matrix(2*n_nodes, k);
	for (int i = 0; i < bandK->m*(2*bandK->k+1); i++) bandK->data[i] = 0.0;
	for (int t = 0; t < n_triplets; t++) bandK->a[triplets[t].i][triplets[t].j] += triplets[t].val;

    double bandRHS[2*n_nodes];
	memcpy(bandRHS, permRHS, 2*n_nodes*sizeof(double));

	// --- 4 ---
	SymBandMatrix *symK = allocate_sym_band_matrix(2*n_nodes, k);
	for (int i = 0; i < symK->m*(symK->k+1); i++) symK->data[i] = 0.0;
	for (int t = 0; t < n_triplets; t++) if (triplets[t].i <= triplets[t].j) symK->a[triplets[t].i][triplets[t].j] += triplets[t].val;

	double symRHS[2*n_nodes];
	memcpy(symRHS, permRHS, 2*n_nodes*sizeof(double));

	// Solve linear system

	// --- 1 ---
	lu(K);
	solve(K, RHS);

	// --- 2 ---
	lu(permK);
	solve(permK, permRHS);

	double AUX[2*n_nodes];
	memcpy(AUX, permRHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) permRHS[i] = AUX[perm[i]]; // re-permute RHS

	// --- 3 ---
	lu_band(bandK);
	solve_band(bandK, bandRHS);

	memcpy(AUX, bandRHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) bandRHS[i] = AUX[perm[i]]; // re-permute RHS

	// --- 4 ---
	cholesky(symK);
	solve_cholesky(symK, symRHS);

	memcpy(AUX, symRHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) symRHS[i] = AUX[perm[i]]; // re-permute RHS

	for (int i = 0; i < 2*n_nodes; i++) {
		if (abs(RHS[i] - permRHS[i]) > 1e-10) {
			printf("Error RHS vs permRHS: %d %lf %lf\n", i, RHS[i], permRHS[i]);
		}
		if (abs(permRHS[i] - bandRHS[i]) > 1e-10) {
			printf("Error permRHS vs bandRHS : %d %lf %lf\n", i, permRHS[i], bandRHS[i]);
		}
		if (abs(bandRHS[i] - symRHS[i]) > 1e-10) {
			printf("Error: %d %lf %lf\n", i, bandRHS[i], symRHS[i]);
		}
	}

	// Visualization in Gmsh
	visualize_in_gmsh(symRHS, gmsh_num, n_nodes);

	// Run the Gmsh GUI; Comment this line if you do not want the Gmsh window to launch
	gmshFltkRun(&ierr);

	// Free stuff
	// --- 1 ---
	free(RHS);
	free(gmsh_num);
	free(coord);
	free(triplets);
	free_matrix(K);
	gmshFinalize(&ierr);

	// --- 2 ---
	free_matrix(permK);

	// --- 3 ---
	free_band_matrix(bandK);

	// --- 4 ---
	free_sym_band_matrix(symK);

	return 0;
}