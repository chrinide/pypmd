#!/usr/bin/env python

import os
import sys
import time
import h5py
import numpy
import ctypes

import logger
import param
import data
import chkfile


# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

# Ultra slow in python just testing, TODO: add opt code
def rho_grad(self,point):
    rho = 0.0
    gun = numpy.zeros(self.nmo)
    gun1 = numpy.zeros((self.nmo,3))
    grad = numpy.zeros(3)
    fun = numpy.zeros((3))
    fun1 = numpy.zeros((3))
    # TODO: avoid natm loop and loop only over total shells
    for ic in range(self.natm):
        xcoor = point - self.coords[ic,:]
        dis2 = numpy.einsum('i,i->', xcoor, xcoor)
        for m in range(self.ngroup[ic]):
            if (dis2 >= self.rcutte[m,ic]): continue
            k = self.nuexp[m,0,ic]
            ori = -self.oexp[k-1]
            dp2 = ori + ori
            aexp = numpy.exp(ori*dis2)
            for jj in range(self.nzexp[m,ic]):
                i = self.nuexp[m,jj,ic]
                itip = self.ityp[i-1]-1
                it = data.nlm[itip,:]
                for j in range(3):
                    n = it[j]
                    x = xcoor[j]
                    if (n == 0):
                        dp2x = dp2*x
                        fun1[j] = dp2x
                        fun[j] = 1.0
                    elif (n == 1):
                        x2 = x*x
                        dp2x2 = dp2*x2
                        fun1[j] = 1.0+dp2x2
                        fun[j] = x
                    elif (n == 2):
                        x2 = x*x
                        dp2x2 = dp2*x2
                        fun1[j] = x*(2.0+dp2x2)
                        fun[j] = x2
                    elif (n == 3):
                        x2 = x*x
                        dp2x2 = dp2*x2
                        fun1[j] = x2*(3.0+dp2x2)
                        fun[j] = x*x2
                    elif (n == 4):
                        x2 = x*x
                        dp2x2 = dp2*x2
                        fun1[j] = x2*x*(4.0+dp2x2)
                        fun[j] = x2*x2
                    elif (n == 5):
                        x2 = x*x
                        dp2x2 = dp2*x2
                        fun1[j] = x2*x2*(5.0+dp2x2)
                        fun[j] = x2*x2*x
                f12 = fun[0]*fun[1]*aexp
                f123 = f12*fun[2]
                fa = fun1[0]*fun[1]*fun[2]*aexp
                fb = fun1[1]*fun[0]*fun[2]*aexp
                fc = fun1[2]*f12
                for j in range(self.nmo):
                    cfj = self.mo_coeff[i-1,j]
                    gun[j] += cfj*f123
                    gun1[j,0] += cfj*fa
                    gun1[j,1] += cfj*fb
                    gun1[j,2] += cfj*fc

    rho = numpy.einsum('i,i->', self.mo_occ, gun*gun)
    grad[0] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,0])
    grad[1] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,1])
    grad[2] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,2])
    grad = grad + grad
    gradmod = numpy.linalg.norm(grad)

    return rho, grad, gradmod

class Becke(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
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
        self.points = None
        self.weights = None
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

        logger.info(self,'* Grid Info')
        logger.info(self,'Tot grids %d', len(self.weights))
        logger.info(self,'')
     
        return self

    def build(self):

        t0 = time.clock()
        logger.TIMER_LEVEL = 3
    
        # 1) Read info
        logger.info(self,'Reading HDF5 file')
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
            self.points = f['becke/coords'].value
            self.weights = f['becke/weights'].value
        logger.info(self,'Finish with HDF5 file')

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()
    
        rho = 0
        for i in range(len(self.weights)):
            tmp = rho_grad(self,self.points[i])[0]
            rho += tmp*self.weights[i]
        logger.info(self,'The value of rho is', rho)

        logger.timer(self,'Becke integration done', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    becke = Becke(name)
    becke.verbose = 4
    becke.kernel()

    from pyscf import lib, dft
    from pyscf.dft import numint

    chkname = 'h2o.chk'
    mol = lib.chkfile.load_mol(chkname)
    mf_mo_coeff = lib.chkfile.load(chkname, 'scf/mo_coeff')
    mf_mo_occ = lib.chkfile.load(chkname, 'scf/mo_occ')
    coords = numpy.reshape(becke.points, (-1,3))
    ao = dft.numint.eval_ao(mol, coords, deriv=0)
    rho = dft.numint.eval_rho2(mol, ao, mf_mo_coeff, mf_mo_occ, xctype='LDA')
    rho = numpy.einsum('i,i->',rho,becke.weights)
    print rho

