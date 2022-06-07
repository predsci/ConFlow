######################################
# Conflow
######################################

FFLAGS=-O3 -march=native -ftree-parallelize-loops=${OMP_NUM_THREADS} -I/opt/psi/gnu/ext_deps/deps/hdf5/include

LDFLAGS=-L/opt/psi/gnu/ext_deps/deps/hdf5/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl

conflow : plmcoef.o plm.o four1.o psi_io_h5.o conflow.o
	gfortran ${FFLAGS} plmcoef.o plm.o four1.o psi_io_h5.o conflow.o ${LDFLAGS} -o conflow
	rm *.o

conflow.o : conflow.f90
	gfortran ${FFLAGS} -c conflow.f90

plmcoef.o : plmcoef.f90
	gfortran ${FFLAGS} -c plmcoef.f90

plm.o : plm.f90
	gfortran ${FFLAGS} -c plm.f90

four1.o : four1.f90
	gfortran ${FFLAGS} -c four1.f90

psi_io_h5.o: psi_io_h5.f90
	gfortran ${FFLAGS} -c psi_io_h5.f90

