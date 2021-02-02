all:
	mpif90 -qopenmp -module $(MODDIR) -I/$(INCDIR) nspcg.f PlumeModels.f si3d_types.f90 si3d_ecomod.f90 si3d_utils.f90 si3d_boundaryconditions.f90 si3d_mixing.f90 si3d_procedures.f90 si3d.f90 -c
	mpif90 -qopenmp -o psi3d *.o -module $(MODDIR) -I/$(INCDIR) -L/$(LIBDIR) -lturbulence_prod -lutil_prod
omp:
	ifort -qopenmp -qopenmp-link=static -module $(MODDIR) -I/$(INCDIR) nspcg.f PlumeModels.f si3d_types.f90 si3d_ecomod.f90 si3d_utils.f90 si3d_boundaryconditions.f90 si3d_mixing.f90 si3d_procedures.f90 si3d.f90 -c
	ifort -qopenmp -qopenmp-link=static -o psi3d *.o -module $(MODDIR) -I/$(INCDIR) -L/$(LIBDIR) -lturbulence_prod -lutil_prod
sec:
	ifort -qopenmp -module $(MODDIR) -I/$(INCDIR) nspcg.f PlumeModels.f si3d_types.f90 si3d_ecomod.f90 si3d_boundaryconditions.f90 si3d_mixing.f90 si3d_procedures.f90 si3d.f90 -c -O2
	ifort -qopenmp -o psi3d *.o -module $(MODDIR) -I/$(INCDIR) -L/$(LIBDIR) -lturbulence_prod -lutil_prod -O2
old:
	ifort -I/opt/mpich2/gnu/include -I/opt/mpich2/gnu/include -L/opt/mpich2/gnu/lib -L/opt/mpich2/gnu/lib -lmpichf90 -Wl,-rpath -Wl,/opt/mpich2/gnu/lib -lmpichf90 -lmpich -lopa -lpthread -lrt -luuid -lrt -o psi3d si3d_types.f90 si3d_utils.f90 si3d_mixing.f90 si3d_boundaryconditions.f90 si3d_procedures.f90 nspcg.f PlumeModels.f si3d.f90 -lm -qopenmp -O2
mpi:
	mpif90 -qopenmp -qopenmp-link=static si3d_types.f90 si3d_utils.f90 si3d_mixing.f90 si3d_boundaryconditions.f90 si3d_procedures.f90 nspcg.f PlumeModels.f si3d.f90 -g -O2
mmm:
	mpif90-vt -c -vt:f90 ifort -auto -qopenmp -qopenmp-link=static si3d_types.f90 si3d_utils.f90 si3d_mixing.f90 si3d_boundaryconditions.f90 si3d_procedures.f90 nspcg.f PlumeModels.f si3d.f90 -g -O2
mpich:
	/opt/mpich2-1.2/bin/mpif90 -qopenmp -qopenmp-link=static si3d_types.f90 si3d_utils.f90 si3d_mixing.f90 si3d_boundaryconditions.f90 si3d_procedures.f90 nspcg.f PlumeModels.f si3d.f90 -g -O2
