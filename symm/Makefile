#! /usr/bin/make

include ../make.inc

obj = libsymm.a
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

OBJECTS = mod_sym.o init_sym.o end_sym.o numbersym.o

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

