//#include "function.h"
#include <math.h>
#include "constant.h"
#include<stdlib.h>

void initialize_coords() {
  int k;
  long x, y;
  long n;
  double size_unit_cell[NDIM];
  double coord[NDIM];
  for (k=0; k < NDIM; k++) {
    size_unit_cell[k] = L[k]/NCELL;
  }
  n=0;
  for (x=0; x < NCELL; x++) {
    coord[0] = (x-0.5)*size_unit_cell[0] - LH[0];
    for (y=0; y < NCELL; y++) {
      coord[1] = (y-0.5)*size_unit_cell[1] - LH[1];
    
      for(k=0; k < NDIM; k++) {
        r[k][n] = coord[k];
      }
      n++;
    }
  }
}