#!/usr/bin/env make

CC = gcc
FC = gfortran
LD = gcc
FDEBUG = -fpic -Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
CDEBUG = -fpic -Wpedantic -g -pg -Wunused -Wall
LFLAGS = -shared
FFLAGS = $(FDEBUG) -O3 -mtune=native -fopenmp -fpic
CFLAGS = $(CDEBUG) -O3 -mtune=native -fpic

all: libfapi.so libcapi.so
 
FOBJECTS = mod_prec.o mod_io.o mod_memory.o mod_param.o mod_futils.o \
mod_mole.o mod_basis.o mod_gto.o mod_fields.o mod_surface.o \
surface.o fields.o lebgrid.o  

COBJECTS = csurf.o

libfapi.so: $(FOBJECTS)
	$(LD) $(LFLAGS) -o libfapi.so $(FOBJECTS) -lgfortran

libcapi.so: $(COBJECTS)
	$(LD) $(LFLAGS) -o libcapi.so $(COBJECTS) -lm 

clean:
	/bin/rm -f *.so *.o *.mod *.pyc

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .f .F .c

.f90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.f.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

csurf.o: csurf.h

surface.o surface.mod mod_fields.o mod_fields.mod fields.o fields.mod : mod_basis.mod
surface.o surface.mod mod_surface.o mod_surface.mod mod_odeint.o mod_odeint.mod fields.o fields.mod : mod_fields.mod
mod_surface.o mod_surface.mod mod_mole.o mod_mole.mod mod_basis.o mod_basis.mod : mod_memory.mod
surface.o surface.mod mod_surface.o mod_surface.mod mod_fields.o mod_fields.mod mod_basis.o mod_basis.mod fields.o fields.mod : mod_mole.mod
surface.o surface.mod fields.o fields.mod : mod_param.mod
surface.o surface.mod mod_surface.o mod_surface.mod mod_param.o mod_param.mod mod_odeint.o mod_odeint.mod mod_mole.o mod_mole.mod mod_memory.o mod_memory.mod mod_gto.o mod_gto.mod mod_futils.o mod_futils.mod mod_fields.o mod_fields.mod mod_basis.o mod_basis.mod fields.o fields.mod : mod_prec.mod
surface.o surface.mod : mod_surface.mod
surface.o surface.mod mod_surface.o mod_surface.mod mod_param.o mod_param.mod mod_odeint.o mod_odeint.mod mod_memory.o mod_memory.mod mod_gto.o mod_gto.mod mod_futils.o mod_futils.mod fields.o fields.mod : mod_io.mod
