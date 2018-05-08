#include <stdlib.h>
#include <stdio.h>
#include <xc.h>

void calc_epsxc_PBE( long long Npoints, double *Rhoe, double *gRhoe2, double *epsxc )
{
  double *eps_x = malloc( Npoints*sizeof(double) );
  double *eps_c = malloc( Npoints*sizeof(double) );

  xc_func_type xc_func;

  // PBE exchange 
  xc_func_init( &xc_func, 101, XC_UNPOLARIZED );
  xc_gga_exc( &xc_func, Npoints, Rhoe, gRhoe2, eps_x );
  xc_func_end( &xc_func );

  // VWN correlation
  // LDA_C_VWN_1 = 28
  // LDA_C_VWN   = 7
  xc_func_init( &xc_func, 130, XC_UNPOLARIZED );
  xc_gga_exc( &xc_func, Npoints, Rhoe, gRhoe2, eps_c );
  xc_func_end( &xc_func );

  long long ip;
  for( ip = 0; ip < Npoints; ip++ ) {
    epsxc[ip] = eps_x[ip] + eps_c[ip];
  }

  free( eps_x ); eps_x = NULL;
  free( eps_c ); eps_c = NULL;

}

