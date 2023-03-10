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
	for (int t=0; t<n_triplets; t++) K->a[triplets[t].i][triplets[t].j] += triplets[t].val;

	// --- 2 ---
	int *perm = (int *) malloc(2*n_nodes*sizeof(int));
	compute_permutation(perm, coord, n_nodes, triplets, n_triplets);

	Matrix *permK = allocate_matrix(2*n_nodes, 2*n_nodes);
	for (int i = 0; i < permK->m; i++) for (int j = 0; j < permK->n; j++) permK->a[i][j] = 0.0;

	int k = 0; // compute bandwidth
	for (int t = 0; t < n_triplets; t++) {
		if (abs(triplets[t].j - triplets[t].i) > k) k = abs(triplets[t].j - triplets[t].i);
		permK->a[triplets[t].i][triplets[t].j] += triplets[t].val;
	}
    double *permRHS = (double *) malloc(2*n_nodes*sizeof(double));
	memcpy(permRHS, RHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) permRHS[perm[i]] = RHS[i]; // permute RHS

	// --- 3 ---
	BandMatrix *bandK = allocate_band_matrix(2*n_nodes, k);
	for (int i = 0; i < bandK->m*(2*bandK->k+1); i++) bandK->data[i] = 0.0;
	for (int t = 0; t < n_triplets; t++) bandK->a[triplets[t].i][triplets[t].j] += triplets[t].val;

    double *bandRHS = (double *) malloc(2*n_nodes*sizeof(double));
	memcpy(bandRHS, RHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) bandRHS[perm[i]] = RHS[i]; // permute RHS

	// Solve linear system
	double cpu_time_used;
	clock_t start, end;
	FILE *fd = fopen("time_results.csv", "a");
	fprintf(fd, "%d,", n_nodes);

	// --- 1 ---
	start = clock();
	lu(K);
	solve(K, RHS);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(fd, "%lf,", cpu_time_used);

	// --- 2 ---
	start = clock();
	lu(permK);
	solve(permK, permRHS);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(fd, "%lf,", cpu_time_used);

	double *cpyRHS = (double *) malloc(2*n_nodes*sizeof(double));
	memcpy(cpyRHS, RHS, 2*n_nodes*sizeof(double));
	memcpy(RHS, permRHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) permRHS[i] = RHS[perm[i]]; // re-permute RHS

	// --- 3 ---
	start = clock();
	lu_band(bandK);
	solve_band(bandK, bandRHS);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(fd, "%lf\n", cpu_time_used);

	memcpy(RHS, bandRHS, 2*n_nodes*sizeof(double));
	for (int i = 0; i < 2*n_nodes; i++) bandRHS[i] = RHS[perm[i]]; // re-permute RHS

	for (int i = 0; i < 2*n_nodes; i++) {
		if (abs(permRHS[i] - bandRHS[i]) > 1e-10) {
			printf("Error: %d %lf %lf", i, permRHS[i], bandRHS[i]);
		}
	}

	// Visualization in Gmsh
	// visualize_in_gmsh(permRHS, gmsh_num, n_nodes);

	// Run the Gmsh GUI; Comment this line if you do not want the Gmsh window to launch
	// gmshFltkRun(&ierr);

	// Free stuff
	// --- 1 ---
	free(RHS);
	free(cpyRHS);
	free(gmsh_num);
	free(coord);
	free(triplets);
	free_matrix(K);
	gmshFinalize(&ierr);

	// --- 2 ---
	free(perm);
	free_matrix(permK);
	free(permRHS);

	// --- 3 ---
	free_band_matrix(bandK);
	free(bandRHS);

	return 0;
}