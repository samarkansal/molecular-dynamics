#include "function.h"
#include "constant.h"

void integrate_leapfrog () {
  int k,n;
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      rv[k][n] += ra[k][n]*deltat;
      r[k][n]  += rv[k][n]*deltat;
    }
  }
}