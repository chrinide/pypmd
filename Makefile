#!/usr/bin/env make

FC = gfortran
CC = gcc
LD = gfortran
FDEBUG = -Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
LFLAGS = -shared
CFLAGS = -O3 -fpic -mtune=native -fopenmp
FFLAGS = -O3 -fpic -mtune=native -fopenmp $(FDEBUG) 

all: libfapi.so 
 
FOBJECTS = mod_prec.o mod_io.o mod_memory.o mod_param.o mod_futils.o \
mod_math.o mod_mole.o mod_basis.o mod_fields.o mod_surface.o mod_slm.o \
mod_gaunt.o mod_gto.o mod_eval.o surface.o lebgrid.o

COBJECTS = misc.o

libfapi.so: $(FOBJECTS) $(COBJECTS)
	$(LD) $(LFLAGS) -o libfapi.so $(FOBJECTS) $(COBJECTS)

clean:
	/bin/rm -f *.so *.o *.mod *.pyc

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .f .F .c

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

.f90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.f.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

