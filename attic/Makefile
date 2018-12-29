#!/usr/bin/env make

CC = gcc
LD = gcc
CDEBUG = -fpic -Wpedantic -g -pg -Wunused -Wall
LFLAGS = -shared
CFLAGS = $(CDEBUG) -O3 -mtune=native -fpic

all: libcleb.so
 
COBJECTS = cleb.o

libcleb.so: $(COBJECTS)
	$(LD) $(LFLAGS) -o libcleb.so $(COBJECTS) -lm 

clean:
	/bin/rm -f *.so *.o *.mod *.pyc

.SUFFIXES:
.SUFFIXES: .o .c

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

