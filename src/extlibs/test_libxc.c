#include <stdlib.h>
#include <stdio.h>

void calc_epsxc_VWN( long long Npoints, double *Rhoe, double *epsxc );
void calc_Vxc_VWN( long long Npoints, double *Rhoe, double *V_xc );

int main()
{
  long long Npoints = 5;

  double* Rhoe = malloc( Npoints*sizeof(double) );
  double* V_xc = malloc( Npoints*sizeof(double) );
  double* epsxc = malloc( Npoints*sizeof(double) );

  Rhoe[0] = 0.1;
  Rhoe[1] = 0.2;
  Rhoe[2] = 0.3;
  Rhoe[3] = 0.4;
  Rhoe[4] = 0.5;

  calc_Vxc_VWN( Npoints, Rhoe, V_xc );
  calc_epsxc_VWN( Npoints, Rhoe, epsxc );
  
  int ip;
  printf("\n");
  for( ip = 0; ip < Npoints; ip++ ) {
    printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], V_xc[ip]);
  }

  free( Rhoe ); Rhoe = NULL;
  free( V_xc ); V_xc = NULL;

  return 0;
}
