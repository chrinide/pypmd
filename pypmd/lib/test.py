#!/usr/bin/env python

import numpy, ctypes
from pypmd import lib

libleb = lib.load_library('libsurf')

npang = 5810
libleb.test(ctypes.c_int(npang)) 
