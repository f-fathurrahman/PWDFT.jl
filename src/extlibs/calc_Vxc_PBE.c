#include <stdlib.h>
#include <stdio.h>
#include <xc.h>

void calc_Vxc_PBE( long long Npoints, double *Rhoe, double *gRhoe2,
                   double *V_xc, double *Vg_xc )
{
  double *vrho_x = malloc( Npoints*sizeof(double) );
  double *vrho_c = malloc( Npoints*sizeof(double) );

  double *vgrho_x = malloc( Npoints*sizeof(double) );
  double *vgrho_c = malloc( Npoints*sizeof(double) );

  long long ip;
  for( ip = 0; ip < Npoints; ip++ ) {
    vrho_x[ip] = 0.0;
    vrho_c[ip] = 0.0;
    vgrho_x[ip] = 0.0;
    vgrho_x[ip] = 0.0;
  }
  
  xc_func_type xc_func;

  // PBE exchange 
  xc_func_init( &xc_func, 101, XC_UNPOLARIZED );
  xc_gga_vxc( &xc_func, Npoints, Rhoe, gRhoe2, vrho_x, vgrho_x );
  xc_func_end( &xc_func );

  // PBE correlation
  xc_func_init( &xc_func, 130, XC_UNPOLARIZED );
  xc_gga_vxc( &xc_func, Npoints, Rhoe, gRhoe2, vrho_c, vgrho_c );
  xc_func_end( &xc_func );

  for( ip = 0; ip < Npoints; ip++ ) {
    V_xc[ip] = vrho_x[ip] + vrho_c[ip];
    Vg_xc[ip] = vgrho_x[ip] + vgrho_c[ip];
  }

  free( vrho_x ); vrho_x = NULL;
  free( vrho_c ); vrho_c = NULL;
  free( vgrho_x ); vgrho_x = NULL;
  free( vgrho_c ); vgrho_c = NULL;
}

