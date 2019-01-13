#!/bin/bash
gfortran main.f90 bzip.o mod_io.o -I../ -lbz2 -o main.x
./main.x n2_cas.den.bz2 > n2_cas.den
