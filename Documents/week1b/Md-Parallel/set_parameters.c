#include "function.h"
#include <math.h>
#include "constant.h"
#include <stdlib.h>

void set_parameters() {
  int k;
  for(k=0; k < NDIM; k++) {
    r[k]  = (double *)malloc(NATOM*sizeof(double));
    rv[k] = (double *)malloc(NATOM*sizeof(double));
    ra[k] = (double *)malloc(NATOM*sizeof(double));
  }
  rcut = pow(2.0, 1.0/6.0);
  for (k=0; k < NDIM; k++) {
    L[k]  = NCELL/(sqrt(density));
    LH[k] = 0.5*L[k];
  }
  VMAG = sqrt(temperature*NDIM*(1.0 - 1.0/NATOM)); 
}
