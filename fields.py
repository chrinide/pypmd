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
import misc


# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str
    

libfapi = misc.load_library('libfapi.so')

    
def density_grad_shell(self,point):
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


def density_grad(self,point):
    rho = 0.0
    grad = numpy.zeros(3)
    gun = numpy.zeros(self.nmo)
    gun1 = numpy.zeros((self.nmo,3))
    fun = numpy.zeros((3))
    fun1 = numpy.zeros((3))
    for i in range(self.nprims):
        ic = self.icen[i]-1
        itip = self.ityp[i]-1
        it = data.nlm[itip,:]
        ori = -self.oexp[i]
        dp2 = ori+ori
        xcoor = point - self.coords[ic,:]
        dis2 = numpy.einsum('i,i->', xcoor, xcoor)
        aexp = numpy.exp(ori*dis2)
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

        f23 = fun[1]*fun[2]*aexp
        f13 = fun[0]*fun[2]*aexp
        f12 = fun[0]*fun[1]*aexp
        for j in range(self.nmo):
            cfj = self.mo_coeff[i,j]
            gun[j] += cfj*fun[0]*f23
            gun1[j,0] += cfj*(fun1[0]*f23)
            gun1[j,1] += cfj*(fun1[1]*f13)
            gun1[j,2] += cfj*(fun1[2]*f12)
         
    rho = numpy.einsum('i,i->', self.mo_occ, gun*gun)
    grad[0] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,0])
    grad[1] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,1])
    grad[2] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,2])
    grad = grad + grad
    gradmod = numpy.linalg.norm(grad)

    return rho, grad, gradmod


class Fields(object):

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
        self.nshells = None
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
        logger.info(self,'Total number of shells %d' % self.nshells)
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
            self.nshells = f['basis/nshells'].value

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        logger.info(self,'Finish with HDF5 file')
        logger.timer(self,'Info readed', t0)
        logger.info(self,'')

        return self

    def density_grad(self,point):
        return density_grad(self,point) 
    def density_grad_shell(self,point):
        return density_grad_shell(self,point) 
    def f_density_grad_test(self,point):
        libfapi.driver(ctypes.c_int(self.nmo),  
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
                       self.coords.ctypes.data_as(ctypes.c_void_p)) 
        return self

    kernel = build

if __name__ == '__main__':
    point = [0,0,0.3]
    name = 'h2o.wfn.h5'
    wfn = Fields(name)
    wfn.verbose = 4
    wfn.kernel()
    valores = wfn.density_grad(point)
    print valores
    valores = wfn.density_grad_shell(point)
    print valores
    wfn.f_density_grad_test(point)

    from pyscf import lib, dft
    from pyscf.dft import numint

    chkname = 'h2o.chk'
    mol = lib.chkfile.load_mol(chkname)
    mf_mo_coeff = lib.chkfile.load(chkname, 'scf/mo_coeff')
    mf_mo_occ = lib.chkfile.load(chkname, 'scf/mo_occ')
    a = point
    a = numpy.asarray(a)
    a = numpy.reshape(a, (-1,3))
    ao = dft.numint.eval_ao(mol, a, deriv=1)
    rho = dft.numint.eval_rho2(mol, ao, mf_mo_coeff, mf_mo_occ, xctype='GGA')
    gradmod = numpy.linalg.norm(rho[-3:])
    print('Rho point = %s %s ' % (rho,gradmod))

