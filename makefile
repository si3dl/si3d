FC=ifort
pFC=mpif90

psi3d: si3d_utils.o si3d_boundaryconditions.o si3d_ecomod.o si3d_mixing.o si3d_procedures.o si3d_types.o si3d.o nspcg.o PlumeModels.o
	$(FC)  $^ -o $@ -O2 -fopenmp

si3d.o:si3d.f90 si3d_procedures.o
	$(FC) -c $< -qopenmp -qopenmp-link=static -o $@ -O2

si3d_procedures.o:si3d_procedures.f90 si3d_types.o si3d_boundaryconditions.o si3d_mixing.o si3d_ecomod.o
	$(FC) -c $< -qopenmp -qopenmp-link=static -o $@ -O2

si3d_mixing.o:si3d_mixing.f90 si3d_types.o #turbulence.o kpp.o mtridiagonal.o eqstate.o

si3d_boundaryconditions.o:si3d_boundaryconditions.f90 si3d_types.o si3d_utils.o
	$(FC) -c $< -qopenmp -qopenmp-link=static -o $@ -O2

si3d_utils.o:si3d_utils.f90 si3d_ecomod.o si3d_types.o
	$(FC) -c $< -qopenmp -qopenmp-link=static -o $@ -O2 

si3d_ecomod.o:si3d_ecomod.f90 si3d_types.o

%.o:%.f90
	$(FC) -c $< -o $@ -O2

%.o:%.f
	$(FC) -c $< -o $@ -O2

%.mod:%.f90
	$(FC) -c $< -o $@ -O2

.PHONY:clean

clean:
	rm *.o
	rm psi3d
	rm *.mod
