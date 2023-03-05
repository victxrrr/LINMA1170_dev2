#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	Matrix * K = allocate_matrix(2*n_nodes, 2*n_nodes);
	for (int t=0; t<n_triplets; t++){
		K->a[triplets[t].i][triplets[t].j] += triplets[t].val;
	}

	// ----------------- Write matrix to file -----------------
	FILE *fd1 = fopen("matrix.csv", "w+");
	for (int i=0; i<K->m; i++){
		for (int j=0; j<K->n; j++){
			if (j==K->n-1) fprintf(fd1, "%f", K->a[i][j]);
			else fprintf(fd1, "%f,", K->a[i][j]);
		}
		fprintf(fd1, "\n");
	}
	fclose(fd1);
	// --------------------------------------------------------

	// Solve linear system
	lu(K);
	solve(K, RHS);

	// ----------------- Write solution to file -----------------
	FILE *fd2 = fopen("solution.csv", "w+");
	for (int i=0; i<K->m; i++){
		for (int j=0; j<K->n; j++){
			if (j==K->n-1) fprintf(fd2, "%f", K->a[i][j]);
			else fprintf(fd2, "%f,", K->a[i][j]);
		}
		fprintf(fd2, "\n");
	}
	fclose(fd2);
	// --------------------------------------------------------

	// Visualization in Gmsh
	visualize_in_gmsh(RHS, gmsh_num, n_nodes);

	// Run the Gmsh GUI; Comment this line if you do not want the Gmsh window to launch
	gmshFltkRun(&ierr);

	// Free stuff
	free(RHS);
	free(gmsh_num);
	free(coord);
	free(triplets);
	free_matrix(K);
	gmshFinalize(&ierr);
	return 0;
}