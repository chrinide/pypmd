#! /usr/bin/make

include make.inc

all: promolden
exe = pmd
obj = mod_prec.o mod_io.o mod_param.o mod_utils.o mod_linalg.o \
			mod_memory.o mod_datatm.o mod_wfn.o mod_fields.o mod_odeint.o \
			mod_surf.o lebgrid.o quadratures.o mod_quad.o promolden.o

.SUFFIXES: .o .f90 .f

%.o: %.F90
	$(FC) -c $(FFLAGS) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.f90
	$(FC) -c $(FFLAGS) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.f
	$(FC) -c $(FFLAGS) $(OPT) $(FDEBUG) -o $@ $<

promolden: $(obj)
	$(FC) $(LDFLAG) -o $(exe) $(obj) 

clean:
	@rm -f $(obj) *.mod $(exe)

.PHONY: all promolden clean

mod_surf.o mod_surf.mod : mod_fields.mod
mod_wfn.o mod_wfn.mod : mod_linalg.mod
mod_wfn.o mod_wfn.mod mod_surf.o mod_surf.mod : mod_memory.mod
promolden.o promolden.mod mod_wfn.o mod_wfn.mod mod_surf.o mod_surf.mod mod_datatm.o mod_datatm.mod : mod_param.mod
promolden.o promolden.mod mod_wfn.o mod_wfn.mod mod_utils.o mod_utils.mod mod_surf.o mod_surf.mod mod_param.o mod_param.mod mod_odeint.o mod_odeint.mod mod_memory.o mod_memory.mod mod_linalg.o mod_linalg.mod mod_fields.o mod_fields.mod mod_datatm.o mod_datatm.mod : mod_prec.mod
promolden.o promolden.mod : mod_surf.mod
promolden.o promolden.mod mod_surf.o mod_surf.mod mod_fields.o mod_fields.mod : mod_wfn.mod
promolden.o promolden.mod mod_wfn.o mod_wfn.mod mod_utils.o mod_utils.mod mod_surf.o mod_surf.mod mod_param.o mod_param.mod mod_memory.o mod_memory.mod mod_linalg.o mod_linalg.mod : mod_io.mod
