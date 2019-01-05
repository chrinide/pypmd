#!/usr/bin/env python

import numpy

MAX_MEMORY = 4000
TMPDIR = '.' 

L_MAX = 5
ANGULAR = 'spdfgh'
ANGULARMAP = {'s': 0,
              'p': 1,
              'd': 2,
              'f': 3,
              'g': 4,
              'h': 5}
        
MGRP = 200
NGTOH = 21
LTYPE = ['S', 'P', 'D', 'F', 'G', 'H']
NCART = [1, 3, 6, 10, 15, 21]

VERBOSE_DEBUG  = 5
VERBOSE_INFO   = 4
VERBOSE_NOTICE = 3
VERBOSE_WARN   = 2
VERBOSE_ERR    = 1
VERBOSE_QUIET  = 0
VERBOSE_CRIT   = -1
VERBOSE_ALERT  = -2
VERBOSE_PANIC  = -3

BOHR = 0.52917721092  # Angstroms

GRADEPS = 1e-10
RHOEPS = 1e-10
MINSTEP = 1e-6
MAXSTEP = 0.75
SAFETY = 0.9
ENLARGE = 1.6

OCCDROP = 1e-6
HMINIMAL = numpy.finfo(numpy.float64).eps

