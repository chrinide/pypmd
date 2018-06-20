#! /usr/bin/make

include ../make.inc

obj = libgeom.a
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

OBJECTS = mod_geom.o polich.o cmcq.o loadgeom.o infogeom.o \
					init_geom.o connect.o

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

