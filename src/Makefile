######################################
# Conflow
######################################
###############

conflow : conflow.o plmcoef.o plm.o four1.o
	gfortran -o conflow conflow.o plmcoef.o plm.o four1.o

conflow.o : conflow.f90
	gfortran -c conflow.f90

plmcoef.o : plmcoef.f90
	gfortran -c plmcoef.f90

plm.o : plm.f90
	gfortran -c plm.f90

four1.o : four1.f90
	gfortran -c four1.f90
