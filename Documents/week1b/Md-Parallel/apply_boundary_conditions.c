#include "function.h"
#include <math.h>
#include "constant.h"
#include<stdlib.h>

void apply_boundary_conditions() {
  int k,n;
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      if (r[k][n] > LH[k]) {
        r[k][n] = r[k][n] - L[k];
      }
      if (r[k][n] < -LH[k]) {
        r[k][n] = r[k][n] + L[k];
      }
    }
  }
}
