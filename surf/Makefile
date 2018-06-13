#! /usr/bin/make

include ../make.inc

obj = libsurf.a
all : $(obj)

.SUFFIXES: .o .f90 .f

%.o: %.f
	$(FCOMPL) -I../common  -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

%.o: %.f90
	$(FCOMPL) -I../common  -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -I../common -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

OBJECTS = weightheta.o

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

