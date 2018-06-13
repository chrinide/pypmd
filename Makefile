#! /usr/bin/make

include make.inc

all: promolden
obj: promolden
files = promolden.o

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -Icommon -c $(FFLAGC) -O2 $(FDEBUG) -o $@ $<

libcommon:
	cd common/ && FC="$(FCOMPL)" FFLAGS="$(FFLAGC)" $(MAKE)

promolden: libcommon $(files)
	$(FCOMPL) $(LDFLAG) -o promolden $(files) common/libcommon.a

clean:
	cd common && make clean
	@rm -f $(files) *.mod $(obj)

.PHONY: all promolden libcommon clean
