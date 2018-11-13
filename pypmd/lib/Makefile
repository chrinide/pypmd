#!/usr/bin/env make

FC = gfortran
LD = gfortran
FDEBUG = -Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
FFLAGS = $(FDEBUG) -pipe -falign-commons -mtune=native -fopenmp -fpic 
LFLAGS = -shared
FOPT   = -O3
#-fstack-arrays -ftree-loop-optimize -ffp-contract=fast \
#-ftree-vectorize -ftree-loop-distribution \
#-Ofast -ffast-math -mtune=native -march=native \
#-funsafe-loop-optimizations -funsafe-math-optimizations

all: libsurf.so

OBJECTS = mod_prec.o mod_memory.o mod_io.o \
mod_param.o mod_wfn.o lebgrid.o wfn.o surf.o

libsurf.so: $(OBJECTS)
	$(LD) $(LFLAGS) -o libsurf.so $(OBJECTS)

clean:
	/bin/rm *.so *.o

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .f .F

.f90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F90.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.f.o:
	$(FC) -c $(FFLAGS) -o $@ $<

.F.o:
	$(FC) -c $(FFLAGS) -o $@ $<

mod_wfn.o mod_wfn.mod : mod_memory.mod
wfn.o wfn.mod surf.o surf.mod mod_fields.o mod_fields.mod : mod_param.mod
wfn.o wfn.mod mod_wfn.o mod_wfn.mod mod_param.o mod_param.mod mod_memory.o mod_memory.mod mod_fields.o mod_fields.mod : mod_prec.mod
wfn.o wfn.mod mod_fields.o mod_fields.mod : mod_wfn.mod
mod_wfn.o mod_wfn.mod mod_param.o mod_param.mod mod_memory.o mod_memory.mod mod_fields.o mod_fields.mod : mod_io.mod

