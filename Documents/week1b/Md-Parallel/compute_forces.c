#include "function.h"
#include <math.h>
#include "constant.h"
#include<stdlib.h>
#include <stdio.h>

void compute_forces() {
  long n, k;
  long n1, n2;
  
  double dr[NDIM], fmag, f, rrcut, rr, rri, rri3;
  rrcut = pow(rcut, 2.0);
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      ra[k][n] = 0.0;
    }
  }
  usum       = 0.0;
  virial_sum = 0.0;
  
 for (n1=0; n1 < NATOM; n1++) {
  for (n2=0; n2 < NATOM; n2++) {
    if (n1 < n2) {
      rr=0.0;
      for (k=0; k < NDIM; k++) {
        dr[k] = r[k][n1] - r[k][n2];
          if (fabs(dr[k]) > LH[k]) {
            if (dr[k] > 0) {
              dr[k] = (dr[k] - L[k]);
            } else {
              dr[k] = L[k] + dr[k];
            }
          }
          rr += dr[k]*dr[k];
        }
        if (rr == 0.0) {
          exit(0);
        }
        if (rr <  rrcut) {
          rri   =  1.0/rr;
          rri3  = rri*rri*rri;
          fmag  = 48.0*(rri3*(rri3-0.5)*rri);
          for (k=0; k < NDIM; k++) {
            f = fmag*dr[k];
            ra[k][n1] += f;
            ra[k][n2] -= f;
          }
          usum       += 4.0*(rri3)*(rri3-1.0) + 1.0;
          virial_sum += fmag*rr;
        }
      }
    }
  }
}