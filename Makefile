#! /usr/bin/make

include make.inc

all: promolden
exe = promolden
obj = promolden.o

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -Icommon -Iwfn -Isurf -c $(LDFLAG) $(FDEBUG) -o $@ $<

libcommon:
	cd common/ && $(MAKE)

libwfn:
	cd wfn/ && $(MAKE)

libsurf:
	cd surf/ && $(MAKE)

promolden: libcommon libwfn libsurf $(obj)
	$(FCOMPL) $(LDFLAG) -o $(exe) $(obj) wfn/libwfn.a surf/libsurf.a \
  common/libcommon.a

clean:
	cd common && make clean
	cd surf && make clean
	cd wfn && make clean
	@rm -f $(obj) *.mod $(exe)

.PHONY: all promolden libcommon clean
