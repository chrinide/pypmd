#! /usr/bin/make

include ../make.inc

obj = libsurf.a
all : $(obj)

.SUFFIXES: .o .f90 .f

%.o: %.f
	$(FCOMPL) -I../common -I../wfn -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.f90
	$(FCOMPL) -I../common -I../wfn -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -I../common -I../wfn -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

OBJECTS = mod_surf.o init_surf.o findnuc.o gradrho.o \
					weightheta.o iscp.o infosurf.o end_surf.o \
					optssurf.o mod_odeint.o ray.o surface.o

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

