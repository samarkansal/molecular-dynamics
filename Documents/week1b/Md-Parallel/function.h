void set_parameters();
//void set_up_job();

void iterate_onestep();
void compute_forces();
void integrate_leapfrog();
void apply_boundary_conditions();
void evaluate_props();
void write2file(long t);