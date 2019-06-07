//#include "function.h"
#include <math.h>
#include "constant.h"
#include<stdlib.h>
#include <stdio.h>

void initialize_vels() {
  long n;
  int k;
  double random;
  double vsum[NDIM];
  for (n=0; n < NATOM; n++) {
    random = drand48();
//     for (k=0; k < NDIM; k++) {
      rv[0][n] = VMAG*cos(2.0*PI*random);
      rv[1][n] = VMAG*sin(2.0*PI*random);
//     }
  }
  for (k=0; k < NDIM; k++) {
    vsum[k] = 0.0;
  }
  
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      vsum[k] += rv[k][n]; 
    }
  }
  
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      rv[k][n] = rv[k][n] - vsum[k]/NATOM;
    }
  }
}
