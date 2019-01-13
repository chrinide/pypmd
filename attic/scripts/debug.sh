#!/bin/bash
# Load the executable file with all symbols
file ../bin/pmd_debug
 
# This settings is necessary for complex breakpoint conditions in fortran...
#set language fortran
set language c
set language auto
 
# Write where I am...
echo In directory:
shell echo pwd
echo \n Files in dir:
shell ls -a
 
# clean the desk: delete all previous breakpoints
d b

# Begin
b esi.f
