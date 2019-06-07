#include "function.h"
#include <math.h>
#include "constant.h"
#include <stdlib.h>
#include<stdio.h>

void write2file(long t) {
  long n;
  FILE *fp;
  char name[10000];
  //printf("%ld**",t);
  sprintf(name, "dataFiles/Atom_position_%ld.dat", t);
  fp = fopen(name, "w");
  for (n=0; n < NATOM; n++) {
    fprintf(fp,"%le %le\n",r[0][n], r[1][n]);
  }
  fclose(fp);
//   sprintf(name, "Measurements.dat", t);
  fp = fopen("Measurements.dat", "a");
  fprintf(fp, "%ld %le %le %le %le \n", t, kin_E, pot_E, tot_E, pressure);
  fclose(fp);
}