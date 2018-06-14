#! /usr/bin/make

include ../make.inc

obj = libwfn.a
all : $(obj)

.SUFFIXES: .o .f90 .f

%.o: %.f
	$(FCOMPL) -I../common -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.f90
	$(FCOMPL) -I../common -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -I../common -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

OBJECTS = mod_wfn.o setupwfn.o rdwfn.o filtergto.o etijcalc.o \
					gtogto.o optswfn.o init_wfn.o loadwfn.o end_wfn.o \
					infowfn.o

clean:
	@rm -f $(OBJECTS) *.mod $(obj)

$(obj): $(OBJECTS) 
	$(AR) $@ $(OBJECTS)
	$(RANLIB) $@

