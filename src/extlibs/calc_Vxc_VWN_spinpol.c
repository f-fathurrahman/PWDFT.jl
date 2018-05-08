#include <stdlib.h>
#include <stdio.h>
#include <xc.h>

void calc_Vxc_VWN_spinpol( long long Npoints, double *Rhoe, double *V_xc )
{
  double *vrho_x = malloc( 2*Npoints*sizeof(double) );
  double *vrho_c = malloc( 2*Npoints*sizeof(double) );

  long long ip;
  for( ip = 0; ip < 2*Npoints; ip++ ) {
    vrho_x[ip] = 0.0;
    vrho_c[ip] = 0.0;
  }

  xc_func_type xc_func;

  // LDA exchange 
  xc_func_init( &xc_func, 1, XC_POLARIZED );
  xc_lda_vxc( &xc_func, Npoints, Rhoe, vrho_x );
  xc_func_end( &xc_func );

  // VWN correlation
  // LDA_C_VWN_1 = 28
  // LDA_C_VWN   = 7
  xc_func_init( &xc_func, 7, XC_POLARIZED );
  xc_lda_vxc( &xc_func, Npoints, Rhoe, vrho_c );
  xc_func_end( &xc_func );

  for( ip = 0; ip < 2*Npoints; ip++ ) {
    V_xc[ip] = vrho_x[ip] + vrho_c[ip];
  }

  free( vrho_x ); vrho_x = NULL;
  free( vrho_c ); vrho_c = NULL;

}

