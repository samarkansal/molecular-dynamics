#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "function.h"
#include "constant.h"
#include "initialize_coords.h"
#include "initialize_vels.h"
/*
#define ntimesteps 100000
#define time_interval 10
#define NDIM 2
#define NATOM 400
#define PI M_PI
#define density 0.8
#define NCELL 20
// #define VMAG 1.0
#define temperature 1.0
#define deltat 0.005

double *r[NDIM], *rv[NDIM], *ra[NDIM];

double rcut;
double L[NDIM];
double LH[NDIM];
double VMAG;
double virial_sum;
double usum;

double kin_E, pot_E, tot_E, pressure;
*/
#define noth 32

#include<omp.h>

void main() {
  long t;
  set_parameters();
  //set_up_job();
  initialize_coords();
  initialize_vels();
  double run_time ;

  double start_time = omp_get_wtime();
  omp_set_num_threads(noth);
#pragma omp parallel num_threads(noth)
  {
  long id=omp_get_thread_num();  
  long t;
  for (t=id; t < ntimesteps; t=t+noth) {
    printf("t=%ld by thread:%ld\n",t,id);
    iterate_onestep();
    evaluate_props();
    if ((t%time_interval) == 0) {
     // printf("+++\n");
      write2file(t);
    }
  }
  
}
run_time = omp_get_wtime() - start_time;
  printf("run time=%f",run_time);
}