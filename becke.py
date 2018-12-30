#!/usr/bin/env python

import os
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
import grids
import misc

signal.signal(signal.SIGINT, signal.SIG_DFL)

libfapi = misc.load_library('libfapi')
libgapi = misc.load_library('libgapi')

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

class Becke(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
        self.small_rho_cutoff = 1e-7
        self.grids = grids.Grids(self.chkfile)
        self.gpu = False
##################################################
# don't modify the following attributes, they are not input options
        self.natm = None
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
            logger.info(self,'Nuclei %d position : %.6f  %.6f  %.6f', i, *self.coords[i])

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

        return self

    def build(self):

        t0 = time.clock()
        logger.TIMER_LEVEL = 3
    
        # 1) Build grid
        self.grids.verbose = self.verbose
        self.grids.stdout = self.stdout
        self.grids.max_memory = self.max_memory
        self.grids.scratch = self.scratch
        self.grids.build()

        # 2) Read info
        with h5py.File(self.chkfile) as f:
            self.natm = f['molecule/natm'].value
            self.coords = f['molecule/coords'].value
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

        npoints = len(self.grids.weights)
        if (self.gpu):
            feval = 'gpu_info'
            drv = getattr(libgapi, feval)
            drv()
        else:
            feval = 'becke_driver'
            drv = getattr(libfapi, feval)
            drv(ctypes.c_int(self.nmo), 
                ctypes.c_int(self.nprims),  
                self.icen.ctypes.data_as(ctypes.c_void_p), 
                self.ityp.ctypes.data_as(ctypes.c_void_p), 
                self.oexp.ctypes.data_as(ctypes.c_void_p), 
                self.ngroup.ctypes.data_as(ctypes.c_void_p), 
                self.nzexp.ctypes.data_as(ctypes.c_void_p), 
                self.nuexp.ctypes.data_as(ctypes.c_void_p), 
                self.rcutte.ctypes.data_as(ctypes.c_void_p), 
                self.mo_coeff.ctypes.data_as(ctypes.c_void_p), 
                self.mo_occ.ctypes.data_as(ctypes.c_void_p), 
                ctypes.c_int(self.natm),  
                self.coords.ctypes.data_as(ctypes.c_void_p), 
                ctypes.c_int(npoints),  
                self.grids.coords.ctypes.data_as(ctypes.c_void_p),
                self.grids.weights.ctypes.data_as(ctypes.c_void_p))
        
        logger.timer(self,'Becke integration done', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    becke = Becke(name)
    becke.verbose = 4
    becke.grids.level = 3
    becke.kernel()

    from pyscf import lib, dft
    from pyscf.dft import numint

    chkname = 'h2o.chk'
    mol = lib.chkfile.load_mol(chkname)
    mf_mo_coeff = lib.chkfile.load(chkname, 'scf/mo_coeff')
    mf_mo_occ = lib.chkfile.load(chkname, 'scf/mo_occ')
    coords = numpy.reshape(becke.grids.coords, (-1,3))
    ao = dft.numint.eval_ao(mol, coords, deriv=0)
    rho = dft.numint.eval_rho2(mol, ao, mf_mo_coeff, mf_mo_occ, xctype='LDA')
    rho = numpy.einsum('i,i->',rho,becke.grids.weights)
    print rho

