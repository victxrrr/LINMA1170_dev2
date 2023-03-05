#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../gmsh-sdk/include/gmshc.h"
#include <math.h>
#include "elasticity.h"


static void p1_stiffness_matrix (double dphi[6], double h[3][3] , double det, double S[6][6]){

  double d1dx = dphi[0];
  double d1dy = dphi[1];
  double d2dx = dphi[2];
  double d2dy = dphi[3];
  double d3dx = dphi[4];
  double d3dy = dphi[5];
  double BT[3][6] = {{d1dx, d2dx, d3dx, 0, 0, 0},
		     {0, 0, 0, d1dy, d2dy, d3dy},
		     {d1dy, d2dy, d3dy, d1dx, d2dx, d3dx}};
  double Bh[6][3];
  for (int i=0;i < 6; i++){
    for (int j=0;j < 3; j++){
      Bh[i][j] = 0;
      for (int k=0;k < 3; k++){
	Bh[i][j] += BT[k][i]*h[k][j];
      }
    }
  }

  for (int i=0;i < 6; i++){
    for (int j=0;j < 6; j++){
      S[i][j] = 0;
      for (int k=0;k < 3; k++){
	S[i][j] += det*0.5*Bh[i][k]*BT[k][j];
      }
    }
  }
}

static void hooke_plain_stress (double E, double nu, double h[3][3]){
  double C = E/(1-nu*nu);
  h[0][0] = h[1][1] = C;
  h[2][2] = C*(1-nu);
  h[0][1] = h[1][0] = C*nu;
  h[0][2] = h[2][0] = h[1][2] = h[2][1] = 0;
}

void p1_stiffness_matrix_plane_stress (double E, double nu, double dphi[6], double det, double S[6][6]){
  double HookeMatrix[3][3];
  hooke_plain_stress (E,nu,HookeMatrix);    
  p1_stiffness_matrix (dphi,HookeMatrix,det,S);  
}

void p1_geometry(const double *x, double *detptr, double dxidx[2][2], double *dphi)
{
  double dxdxi[2][2] = {
    {x[1*2+0]-x[0*2+0],x[2*2+0]-x[0*2+0]},
    {x[1*2+1]-x[0*2+1],x[2*2+1]-x[0*2+1]}
  };
  double det = dxdxi[0][0]*dxdxi[1][1]-dxdxi[0][1]*dxdxi[1][0];
  dxidx[0][0] = dxdxi[1][1]/det;
  dxidx[0][1] = -dxdxi[0][1]/det;
  dxidx[1][0] = -dxdxi[1][0]/det;
  dxidx[1][1] = dxdxi[0][0]/det;
  for (int i = 0; i < 2; ++i) {
    dphi[0*2+i] = -dxidx[0][i] - dxidx[1][i];
    dphi[1*2+i] = dxidx[0][i];
    dphi[2*2+i] = dxidx[1][i];
  }
  *detptr = det;
}

int get_triangles_p1 (size_t ** elementTags, size_t * elementTags_n,
		       size_t ** nodeTags, size_t * nodeTags_n) {
  int ierr;
  gmshModelMeshGetElementsByType(2, elementTags, elementTags_n,
				 nodeTags, nodeTags_n, -1, 0, 1, &ierr);
  return ierr;
}

int get_nodes (size_t ** nodeTags, size_t * nodeTags_n,
		double ** coord, size_t * coord_n) {
  int ierr;
  gmshModelMeshGetNodes(nodeTags, nodeTags_n, coord, coord_n, 0, 0, 2, -1, 1, 0, &ierr);
  return ierr;
}

Matrix * compute_stiffness (int nNodes,
				       double *coord,
				       size_t ntriangles,
				       size_t *triangleNodes){
  double E = 210e9;
  double nu = 0.3;
  Matrix *K = allocate_matrix(2*nNodes, 2*nNodes);
  for (size_t i = 0 ; i < ntriangles ; i++){
    size_t *n = & triangleNodes [3*i];
    double x[6] = {coord[2*(n[0]-1)],coord[2*(n[0]-1)+1],
		   coord[2*(n[1]-1)],coord[2*(n[1]-1)+1],
		   coord[2*(n[2]-1)],coord[2*(n[2]-1)+1]};
    
    double det,dxidx[2][2],dphi[6];
    p1_geometry(x,&det, dxidx, dphi);
    double StiffnessMatrix[6][6];
    p1_stiffness_matrix_plane_stress (E,nu,dphi,det,StiffnessMatrix);
    
    size_t dofs[6] = {2*(n[0]-1),2*(n[1]-1),2*(n[2]-1),
		      2*(n[0]-1)+1,2*(n[1]-1)+1,2*(n[2]-1)+1};
    for (size_t j = 0;j<6;j++){
      for (size_t k = 0;k<6;k++){
        K->a[dofs[j]][dofs[k]] += StiffnessMatrix[j][k];
	// add_to_matrix (K, dofs[j], dofs[k], StiffnessMatrix[j][k]);
      }
    }
  }
  return K;
}

void boundary_conditions (Matrix *K, double *RHS){
  // boundary conditions -- physical names allowed are
  int *dimTags = NULL, ierr;
  size_t dimTags_n;
  gmshModelGetPhysicalGroups(&dimTags, &dimTags_n, -1, &ierr);
  for (size_t i=0;i<dimTags_n;i+=2){  
    char *name =malloc(525*sizeof(char));
    gmshModelGetPhysicalName(dimTags[i],dimTags[i+1],&name, &ierr);
    size_t * nTags = NULL, nTags_n, cc_n;
    double *cc = NULL;
    gmshModelMeshGetNodesForPhysicalGroup(dimTags[i],dimTags[i+1],
					  &nTags, &nTags_n, &cc, &cc_n, &ierr);
    // "clamped" (force x and y = 0)
    if (strcmp(name,"clamped") == 0){
      for (size_t j=0;j<nTags_n;j++){
	size_t r = 2*(nTags[j]-1);
	for (size_t k=0;k<K->n;k++)
    K->a[k][r] = K->a[r][k] = 0;
  K->a[r][r] = 1;
	r = 2*(nTags[j]-1)+1;
	for (size_t k=0;k<K->n;k++)
    K->a[k][r] = K->a[r][k] = 0;
  K->a[r][r] = 1;
      }      
      // printf("clamped found physical group %d\n",dimTags[i+1]);
    }

    // "forcex" (unit force along x)
    if (strcmp(name,"forcex") == 0){
      for (size_t j=0;j<nTags_n;j++){
	size_t r = 2*(nTags[j]-1);
	RHS[r] += 100000;
      }
      // printf("force x found physical group %d\n",dimTags[i+1]);
    }

    // "forcey" (unit force along x)
    if (strcmp(name,"forcey") == 0){
      for (size_t j=0;j<nTags_n;j++){
	size_t r = 2*(nTags[j]-1)+1;
	RHS[r] += 100000;
      }
      // printf("force y found physical group %d\n",dimTags[i+1]);
    }
    if (cc) free(cc);
    if (nTags) free(nTags);                
    free (name);
  }
  if(dimTags)free(dimTags);
}

void assemble_system(Triplet** triplets, int* n_triplets, double** RHS, double** coord, size_t** gmsh_num, int* n_nodes){
  int ierr;
  gmshModelMeshRebuildNodeCache(1,&ierr);
  size_t *triangleTags=0, ntriangles, *triangleNodes=NULL, temp;
  get_triangles_p1 (&triangleTags, &ntriangles,&triangleNodes, &temp);

  size_t *nodeTags=0, nNodes;
  double *coord_bad = 0;
  get_nodes (&nodeTags, &nNodes, &coord_bad, &temp);
  
  *gmsh_num = malloc(nNodes*sizeof(size_t));
  
  *coord = (double*)malloc(2*nNodes*sizeof(double));
  for (size_t i=0;i<nNodes;i++){
    (*coord)[2*(nodeTags[i]-1)] = coord_bad[3*i];
    (*coord)[2*(nodeTags[i]-1)+1] = coord_bad[3*i+1];
    // (*coord)[3*(nodeTags[i]-1)+2] = coord_bad[3*i+2];
    (*gmsh_num)[nodeTags[i]-1] = nodeTags[i];
  }
  Matrix *K = compute_stiffness (nNodes,*coord,ntriangles,triangleNodes);
  *RHS = (double*) malloc (2*nNodes*sizeof(double));
  for (size_t i=0;i<2*nNodes;i++)(*RHS)[i] = 0;
  boundary_conditions (K,*RHS);
  int count = 0; 
  for (size_t i=0; i<K->n; i++){
    for (size_t j=0; j<K->n; j++){
      if (fabs(K->a[i][j]) > 1e-10) count ++;}}

  *triplets = malloc(count*sizeof(Triplet));
  size_t cc = 0;
  for (size_t i=0; i<K->n; i++){
    for (size_t j=0; j<K->n; j++){
      if (fabs(K->a[i][j]) > 1e-10){
        double val = K->a[i][j];
        (*triplets)[cc].i = i;
        (*triplets)[cc].j = j;
        (*triplets)[cc].val = val;
        cc++;
      }
    }
  }


  int n = (int)nNodes;
  *n_nodes = n;
  *n_triplets = count;

  free (coord_bad);
  free (triangleTags);
  free (triangleNodes);
  free_matrix(K);
}

void visualize_in_gmsh(double* SOL, size_t* gmsh_num, int n_nodes){
  int ierr;
  double* sol_3D = malloc(3*n_nodes*sizeof(double)); // in gmsh, vector are 3D
  for (size_t i=0; i<n_nodes; i++){
    sol_3D[3*i+0] = SOL[2*i+0];
    sol_3D[3*i+1] = SOL[2*i+1];
    sol_3D[3*i+2] = 0.;
  }
  char *model_name;
  gmshModelGetCurrent(&model_name, &ierr);
  int view_tag = gmshViewAdd("displacement", -1, &ierr);
  gmshViewAddHomogeneousModelData(view_tag, 0, model_name,"NodeData", gmsh_num, n_nodes, sol_3D, 3*n_nodes, 0., 3, -1, &ierr);
  gmshOptionSetNumber("View.VectorType", 5, &ierr);
  gmshOptionSetNumber("View.DisplacementFactor", 10000, &ierr);
}



