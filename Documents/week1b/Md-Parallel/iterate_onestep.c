#include "function.h"
#include <math.h>
#include "constant.h"
#include<stdlib.h>

void iterate_onestep() {
  compute_forces();
  integrate_leapfrog();
  apply_boundary_conditions();
}