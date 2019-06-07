cc=gcc
objects = main.o set_parameters.o initialize_vels.o iterate_onestep.o compute_forces.o integrate_leapfrog.o apply_boundary_conditions.o evaluate_props.o write2file.o initialize_coords.o
forun : $(objects)
		$(cc) -o forun -fopenmp $(objects) -lm

main.o : main.c function.h constant.h initialize_coords.h initialize_vels.h
	$(cc) -c  main.c -fopenmp
set_parameters.o : set_parameters.c function.h constant.h
	$(cc) -c  set_parameters.c

initialize_coords.o : initialize_coords.c initialize_coords.h constant.h 
	$(cc) -c  initialize_coords.c 
initialize_vels.o : initialize_vels.c  initialize_vels.h constant.h 
	$(cc) -c  initialize_vels.c
iterate_onestep.o : iterate_onestep.c function.h constant.h 
	$(cc) -c  iterate_onestep.c
compute_forces.o : compute_forces.c function.h constant.h 
	$(cc) -c compute_forces.c 
integrate_leapfrog.o : integrate_leapfrog.c function.h constant.h 
	$(cc) -c  integrate_leapfrog.c
apply_boundary_conditions.o : apply_boundary_conditions.c function.h constant.h 
	$(cc) -c  apply_boundary_conditions.c
evaluate_props.o : evaluate_props.c function.h constant.h 
	$(cc) -c  evaluate_props.c
write2file.o : write2file.c function.h constant.h 
	$(cc) -c  write2file.c
clean :
		rm forun $(objects)
