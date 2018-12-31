#!/usr/bin/env python

import sys
import time
import h5py
import numpy
import ctypes
import signal

import logger
import param
import data
import chkfile
import misc
import gto
import tools

signal.signal(signal.SIGINT, signal.SIG_DFL)

libfapi = misc.load_library('libfapi')

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

class Integrals(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
##################################################
# don't modify the following attributes, they are not input options
        self.natm = None
        self.charges = None
        self.symbols = None
        self.coords = None
        self.nmo = None
        self.nprims = None
        self.icen = None
        self.ityp = None
        self.oexp = None
        self.mo_coeff = None
        self.mo_occ = None
        self.ngroup = None 
        self.nzexp = None 
        self.nuexp = None
        self.rcutte = None
        self._keys = set(self.__dict__.keys())

    def dump_input(self):

        if self.verbose < logger.INFO:
            return self

        logger.info(self,'')
        logger.info(self,'******** %s flags ********', self.__class__)
        logger.info(self,'* General Info')
        logger.info(self,'Verbose level %d' % self.verbose)
        logger.info(self,'Scratch dir %s' % self.scratch)
        logger.info(self,'Input h5 data file %s' % self.chkfile)
        logger.info(self,'Max_memory %d MB' % self.max_memory)

        logger.info(self,'* Molecular Info')
        logger.info(self,'Num atoms %d' % self.natm)
        logger.info(self,'Atom Coordinates (Bohr)')
        for i in range(self.natm):
            logger.info(self,'Nuclei %d %s with charge %d position : %.6f  %.6f  %.6f', i, 
                        self.symbols[i], self.charges[i], *self.coords[i])

        logger.info(self,'* Basis Info')
        logger.info(self,'Number of Orbitals %d' % self.nmo)
        logger.info(self,'Number of primitives %d' % self.nprims)
        logger.debug(self,'Number of shells per center %s' % self.ngroup)
        for ic in range(self.natm):
            logger.debug(self,'Basis for center %d' % ic)
            ntfu = numpy.zeros(6, dtype=numpy.int32)
            for m in range(self.ngroup[ic]):
                zz = self.oexp[self.nuexp[m, 0, ic] - 1]
                ii = self.nuexp[m, 0, ic]
                itip = self.ityp[ii - 1] - 1
                it1 = data.nlm[itip, 0]
                it2 = data.nlm[itip, 1]
                it3 = data.nlm[itip, 2]
                isu = it1 + it2 + it3
                ntfu[isu] += 1
                nzicm = self.nzexp[m, ic]
                x1 = self.rcutte[m, ic]
                logger.debug(self, 'Shell Type %s exp %f zero at %f idx %s' % \
                (param.LTYPE[isu],zz,numpy.sqrt(x1),self.nuexp[m, :nzicm, ic]))
        logger.debug(self,'Number of shells of each type %s' % ntfu)
        logger.debug(self,'Ocupation of molecular orbitals %s' % self.mo_occ)
        logger.info(self,'')

        return self

    def build(self):

        t0 = time.clock()
        logger.TIMER_LEVEL = 3
    
        # 1) Read info
        with h5py.File(self.chkfile) as f:
            self.natm = f['molecule/natm'].value
            self.coords = f['molecule/coords'].value
            self.charges = f['molecule/charges'].value
            self.symbols = f['molecule/symbols'].value
            self.nmo = f['basis/nmo'].value
            self.nprims = f['basis/nprims'].value
            self.icen = f['basis/icen'].value
            self.ityp = f['basis/ityp'].value
            self.oexp = f['basis/oexp'].value
            self.mo_coeff = f['basis/mo_coeff'].value
            self.mo_occ = f['basis/mo_occ'].value
            self.ngroup = f['basis/ngroup'].value
            self.nzexp = f['basis/nzexp'].value
            self.nuexp = f['basis/nuexp'].value
            self.rcutte = f['basis/rcutte'].value

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        # 2) Test
        logger.info(self,'Testing orthogonality of orbitals')
        s = gto.eval_overlap(self)
        s = numpy.dot(s,self.mo_coeff)
        s = numpy.dot(self.mo_coeff.T,s)
        tools.dump_tri(self.stdout, s) 
        logger.info(self,'')

        logger.timer(self,'Integrals done', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    integ = Integrals(name)
    integ.verbose = 4
    integ.kernel()

