#! /usr/bin/make

include make.inc

all: promolden
exe = promolden
obj = promolden.o

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -Icommon -Iwfn -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

libcommon:
	cd common/ && FC="$(FCOMPL)" FFLAGS="$(FFLAGC)" $(MAKE)

libwfn:
	cd wfn/ && FC="$(FCOMPL)" FFLAGS="$(FFLAGC)" $(MAKE)

promolden: libcommon libwfn $(obj)
	$(FCOMPL) $(LDFLAG) -o $(exe) $(obj) wfn/libwfn.a surf/libsurf.a \
  common/libcommon.a

clean:
	cd common && make clean
	@rm -f $(obj) *.mod $(exe)

.PHONY: all promolden libcommon clean
