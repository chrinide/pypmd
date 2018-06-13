#! /usr/bin/make

include ../make.inc

obj = libcommon.a
all : $(obj)

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

OBJECTS = mod_prec.o mod_io.o mod_parallel.o mod_param.o 

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

