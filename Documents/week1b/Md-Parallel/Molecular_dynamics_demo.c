#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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

void set_parameters();
void set_up_job();
void initialize_coords();
void initialize_vels();
void iterate_onestep();
void compute_forces();
void integrate_leapfrog();
void apply_boundary_conditions();
void evaluate_props();
void write2file(long t);

void main() {
  long t;
  set_parameters();
  set_up_job();
  
  for (t=0; t < ntimesteps; t++) {
    printf("t=%ld\n",t);
    iterate_onestep();
    evaluate_props();
    if ((t%time_interval) == 0) {
      write2file(t);
    }
  }
}
void write2file(long t) {
  long n;
  FILE *fp;
  char name[10000];
  sprintf(name, "Atom_position_%ld.dat", t);
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
void set_up_job() {
  initialize_coords();
  initialize_vels();
}

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

void iterate_onestep() {
  compute_forces();
  integrate_leapfrog();
  apply_boundary_conditions();
}

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

void integrate_leapfrog () {
  int k,n;
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      rv[k][n] += ra[k][n]*deltat;
      r[k][n]  += rv[k][n]*deltat;
    }
  }
}

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
