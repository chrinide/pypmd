FC = gfortran
FDEBUG = -Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
FFLAGS = -m64 -pipe -falign-commons -fopenmp
CC = gcc
CFLAGS = -O3 -m64 -mtune=native -Wunused
OPT    = 
#-O3 -ftree-loop-optimize \
#-ftree-vectorize -ftree-loop-distribution \
#-Ofast -ffast-math -mtune=native -march=native \
#-funsafe-loop-optimizations -funsafe-math-optimizations
SYSDEP = 
LDFLAG = -m64 -fopenmp
AR = ar rc
RANLIB = ranlib
