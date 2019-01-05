#!/usr/bin/env make

FC = gfortran
LD = gfortran
FDEBUG = -Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
LFLAGS = -shared
FFLAGS = -O3 -fpic -mtune=native -fopenmp $(FDEBUG) 

all: libfapi.so 
 
FOBJECTS = mod_prec.o mod_io.o mod_memory.o mod_param.o mod_futils.o \
mod_math.o mod_mole.o mod_basis.o mod_fields.o mod_surface.o \
surface.o fields.o lebgrid.o

libfapi.so: $(FOBJECTS)
	$(LD) $(LFLAGS) -o libfapi.so $(FOBJECTS)

clean:
	/bin/rm -f *.so *.o *.mod *.pyc

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .f .F

.f90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.f.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

