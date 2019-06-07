#include "function.h"
#include <math.h>
#include "constant.h"
#include<stdlib.h>

void evaluate_props() {
  long k, n;
  double sum_v[NDIM], sum_vv=0.0;
  double v[NDIM];
  
  for (k=0; k < NDIM; k++) {
    sum_v[k] = 0.0;
  }
  
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      v[k]       = rv[k][n] -0.5*deltat*ra[k][n];
      sum_v[k]  += v[k];
      sum_vv += v[k]*v[k];
    }
  }
  kin_E    = 0.5*sum_vv;
  pot_E    = usum;
  tot_E    = kin_E + pot_E;
  pressure = density*(sum_vv + virial_sum)/(NDIM*NATOM); 
}
