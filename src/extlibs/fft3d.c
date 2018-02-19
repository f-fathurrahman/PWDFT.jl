#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

void fftw_fw_fft3d( fftw_complex *in, fftw_complex *out, int nx, int ny, int nz )
{
  fftw_plan plan_forward;
  plan_forward = fftw_plan_dft_3d( nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
  fftw_execute(plan_forward);
  fftw_destroy_plan(plan_forward);
}

void fftw_fw_fft3d_r2c( double *in, fftw_complex *out, int nx, int ny, int nz )
{
  fftw_plan plan_forward;
  plan_forward = fftw_plan_dft_r2c_3d( nx, ny, nz, in, out, FFTW_FORWARD );
  fftw_execute(plan_forward);
  fftw_destroy_plan(plan_forward);
}

void fftw_inv_fft3d( fftw_complex *in, fftw_complex *out, int nx, int ny, int nz )
{
  fftw_plan plan_backward;
  plan_backward = fftw_plan_dft_3d( nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
  fftw_execute(plan_backward);
  fftw_destroy_plan(plan_backward);
  int Ndata = nx*ny*nz;
  double scale = 1.0f/(double)Ndata;
  int inc = 1;
  zdscal_( &Ndata, &scale, out, &inc );
}

void fftw_inv_fft3d_c2r( fftw_complex *in, double *out, int nx, int ny, int nz )
{
  fftw_plan plan_backward;
  plan_backward = fftw_plan_dft_c2r_3d( nx, ny, nz, in, out, FFTW_BACKWARD );
  fftw_execute(plan_backward);
  fftw_destroy_plan(plan_backward);
  int Ndata = nx*ny*nz;
  double scale = 1.0f/(double)Ndata;
  int inc = 1;
  dscal_( &Ndata, &scale, out, &inc );
}
