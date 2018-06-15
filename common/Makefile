#! /usr/bin/make

include ../make.inc

obj = libcommon.a
all : $(obj)

.SUFFIXES: .o .f90 .f

%.o: %.f
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

OBJECTS = mod_prec.o mod_io.o mod_memory.o mod_param.o \
mod_futils.o mod_linalg.o init_param.o quadratures.o \
optsparam.o lebgrid.o

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

