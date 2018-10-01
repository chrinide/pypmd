#FC = gfortran
#FDEBUG = 
##-Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
#FFLAGS = -fopenmp -m64 -pipe -falign-commons 
#CC = gcc
#CFLAGS = -O3 -m64 -mtune=native -Wunused
#OPT    = \
#-O3 -ftree-loop-optimize \
#-ftree-vectorize -ftree-loop-distribution \
#-Ofast -ffast-math -mtune=native -march=native \
#-funsafe-loop-optimizations -funsafe-math-optimizations
#SYSDEP = 
#LDFLAG = -m64 -fopenmp
#AR = ar rc
#RANLIB = ranlib

FC = ifort
FDEBUG = 
FFLAGS = -qopenmp -m64 -static-intel
CC = gcc
CFLAGS = -O3 -m64 
OPT    = -xHOST -O3 -ip -ftz -no-prec-div \
         -no-prec-sqrt -fast-transcendentals
SYSDEP = 
LDFLAG = -m64 -qopenmp
AR = ar rc
RANLIB = ranlib
