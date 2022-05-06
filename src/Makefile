######################################
# Conflow
######################################
###############

FFLAGS=-O3 -march=native -ftree-parallelize-loops=6

conflow : conflow.o plmcoef.o plm.o four1.o
	gfortran ${FFLAGS} -o conflow conflow.o plmcoef.o plm.o four1.o
	rm *.o

conflow.o : conflow.f90
	gfortran ${FFLAGS} -c conflow.f90

plmcoef.o : plmcoef.f90
	gfortran ${FFLAGS} -c plmcoef.f90

plm.o : plm.f90
	gfortran ${FFLAGS} -c plm.f90

four1.o : four1.f90
	gfortran ${FFLAGS} -c four1.f90
