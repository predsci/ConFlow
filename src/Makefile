FC = gfortran

FFLAGS = -O3 -march=native -Wall -Wextra -I/usr/include/hdf5/serial

OBJS = number_types.o plmcoef.o plm.o four1.o psi_io_h5.o zm_parse.o conflow.o

LDFLAGS = -L/usr/lib/x86_64-linux-gnu -lhdf5_serial_fortran -lhdf5_serialhl_fortran -lhdf5_serial -lhdf5_serial_hl

all: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LDFLAGS) -o conflow
	rm *.mod *.o 2>/dev/null

clean:
	rm -f conflow 2>/dev/null
	rm -f *.mod *.o 2>/dev/null

conflow.o : conflow.f90 number_types.f90
	$(FC) $(FFLAGS) -c conflow.f90

plmcoef.o : plmcoef.f90 number_types.f90
	$(FC) $(FFLAGS) -c plmcoef.f90

plm.o : plm.f90 number_types.f90
	$(FC) $(FFLAGS) -c plm.f90

four1.o : four1.f90 number_types.f90
	$(FC) $(FFLAGS) -c four1.f90

psi_io_h5.o: psi_io_h5.f90 number_types.f90
	$(FC) $(FFLAGS) -c psi_io_h5.f90

zm_parse.o: zm_parse.f
	$(FC) $(FFLAGS) -c zm_parse.f

number_types.o: number_types.f90
	$(FC) $(FFLAGS) -c number_types.f90

