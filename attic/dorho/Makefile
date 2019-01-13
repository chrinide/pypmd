#! /usr/bin/make

obj = dorho
all : $(obj)

#FCOMPL = ifort
#FDEBUG = 
#FFLAGC = -m64 -openmp
#CCOMPL = gcc
#CFLAGC = -O2 -m64
#   opt = -xHOST -O3 -ip -unroll -inline -finline-functions \
#         -falign-functions=16 -align all,array16byte \
#				 -no-prec-div -no-prec-sqrt -fast-transcendentals
#SYSDEP =        
#ldflag = $(opt) -m64 -openmp -static-intel  

FCOMPL = gfortran
FDEBUG = -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
FFLAGC = -m64 -fopenmp -pipe -falign-commons
CCOMPL = gcc
CFLAGC = -O2 -m64 -mtune=native -Wunused
opt    = -O3 -funroll-loops -finline -ftree-loop-optimize \
				 -finline-functions -ftree-vectorize -ftree-loop-distribution \
   		   -Ofast -ffast-math -mtune=native -march=native \
         -funsafe-loop-optimizations -funsafe-math-optimizations
SYSDEP = 
ldflag = $(opt) -m64 -fopenmp -static-libstdc++ -static-libgfortran -static-libgcc -lbz2

#FCOMPL = pgfortran -D__PGI
#FDEBUG = -Minfo=accel,inline,ipa,loop,lre,mp,opt,par,unified,vect
#FFLAGC = -m64 -mp  
#CCOMPL = gcc
#CFLAGC = -O2 -m64 -mtune=native -Wunused
#opt    = -fast -fastsse -O4 -Mcache_align
#SYSDEP = 
#ldflag = $(opt) -mp 

# Implicit rules

.SUFFIXES:
.SUFFIXES: .f .o .f90

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) $(opt) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) $(opt) $(FDEBUG) -o $@ $<

%.o: %.f
	$(FCOMPL) -c $(FFLAGC) $(opt) $(FDEBUG) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

# Object files

OBJECTS = mod_prec.o mod_io.o mod_parallel.o mod_futils.o mod_memory.o \
mod_wfn.o mod_param.o mod_bzip.o rdwfn.o setupwfn.o info.o colord.o \
auxcomputerdm.o computerdm.o \
dorho.o

# Clean objects

clean:
	@rm -f core $(OBJECTS) *.mod $(obj)

# Main program target

$(obj): $(SYSDEP) $(VECTOR) $(OBJECTS) 
	$(FCOMPL) $(FFLAGC) -o $(obj) $(OBJECTS) $(SYSDEP) $(ldflag) 

# Include files dependencies

