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