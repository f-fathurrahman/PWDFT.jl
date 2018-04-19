#include <stdio.h>
#include "spglib.h"

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))

int my_spg_find_primitive(double *LatVec, double *atpos, int *types,
                        int Natoms, double symprec)
{
  double lattice[3][3];
  double position[Natoms][3];

  int ia, i;
  for(i = 1; i <= 3; i++) {
    lattice[i-1][0] = LatVec[IDX2F(1,i,3)];
    lattice[i-1][1] = LatVec[IDX2F(2,i,3)];
    lattice[i-1][2] = LatVec[IDX2F(3,i,3)];
    printf("%f %f %f\n", lattice[i-1][0], lattice[i-1][1], lattice[i-1][2]);    
  }

  for(ia = 1; ia <= Natoms; ia++) {
    position[ia-1][0] = atpos[IDX2F(1,ia,3)];
    position[ia-1][1] = atpos[IDX2F(2,ia,3)];
    position[ia-1][2] = atpos[IDX2F(3,ia,3)];
    printf("%3d %3d %f %f %f\n", ia, types[ia-1],
           position[ia-1][0], position[ia-1][1], position[ia-1][2]);
  }

  int num_primitive_atom;
  num_primitive_atom = spg_find_primitive(lattice, position, types, Natoms, symprec);
  
  return num_primitive_atom;
}


void test_2d_array(int a[][3], int dim1)
{
  int i, j;
  for(i = 0; i < dim1; i++) {
    for(j = 0; j < 3; j++) {
      printf("%d ", a[i][j]);
    }
    printf("\n");
  }
  printf("Pass here in test_2d_array\n");
}
