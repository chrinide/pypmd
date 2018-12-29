#!/usr/bin/env python

import sys
import h5py
import time
import numpy
import ctypes
import signal

import data
import misc
import param
import grids
import logger
import chkfile

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

signal.signal(signal.SIGINT, signal.SIG_DFL)

libfapi = misc.load_library('libfapi')
libcapi = misc.load_library('libcapi')


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
    grad = grad/(gradmod+param.HMINIMAL)

    return rho, grad, gradmod


def checkcp(self, x, rho, gradmod):

    iscp = False
    nuc = -2

    for i in range(self.natm):
        r = numpy.linalg.norm(x-self.coords[i])
        if (r < self.epsiscp):
            iscp = True
            nuc = i
            return iscp, nuc

    if (gradmod <= param.GRADEPS):
        iscp = True
        if (rho <= param.RHOEPS): 
            nuc = -1

    return iscp, nuc


def gradrho(self, xpoint, h):

    h0 = h
    niter = 0
    rho, grad, gradmod = rho_grad(self,xpoint)
    grdt = grad
    grdmodule = gradmod

    while (grdmodule > param.GRADEPS and niter < self.mstep):
        niter += 1
        ier = 1
        while (ier != 0):
            xtemp = xpoint + h0*grdt
            rho, grad, gradmod = rho_grad(self,xtemp)
            escalar = numpy.einsum('i,i->',grdt,grad) 
            if (escalar < 0.707):
                if (h0 >= param.MINSTEP):
                    h0 = h0/2.0
                    ier = 1
                else:
                    ier = 0
            else:
                if (escalar > 0.9): 
                    hproo = h0*param.ENLARGE
                    if (hproo < h):
                        h0 = hproo
                    else:
                        h0 = h
                    h0 = numpy.minimum(param.MAXSTEP, h0)
                ier = 0
                xpoint = xtemp
                grdt = grad
                grdmodule = gradmod
            logger.debug(self,'scalar, step in gradrho %.6f %.6f', escalar, h0)

    logger.debug(self,'nsteps in gradrho %d', niter)

    return xpoint, grdmodule


class BaderSurf(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
        self.inuc = 0
        self.epsiscp = 0.180
        self.ntrial = 11
        self.npang = 5810
        self.epsroot = 1e-4
        self.rmaxsurf = 10.0
        self.rprimer = 0.4
        self.backend = 'rkck'
        self.epsilon = 1e-4 
        self.step = 0.1
        self.mstep = 500
        self.csurf = False
##################################################
# don't modify the following attributes, they are not input options
        self.xnuc = None
        self.xyzrho = None
        self.rpru = None
        self.grids = None
        self.rsurf = None
        self.nlimsurf = None
        self.rmin = None
        self.rmax = None
        self.charges = None
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
        logger.info(self,'Input chkfile data file %s' % self.chkfile)
        logger.info(self,'Max_memory %d MB' % self.max_memory)

        logger.info(self,'* Molecular Info')
        logger.info(self,'Num atoms %d' % self.natm)
        logger.info(self,'Atom Coordinates (Bohr)')
        for i in range(self.natm):
            logger.info(self,'Nuclei %d with charge %d position : %.6f  %.6f  %.6f', i, 
                        self.charges[i], *self.coords[i])

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

        logger.info(self,'* Surface Info')
        logger.info(self,'Surface for nuc %d' % (self.inuc+1))
        logger.info(self,'Nuclear coordinate %.6f  %.6f  %.6f', *self.xnuc)
        logger.info(self,'Rmaxsurface %.6f' % self.rmaxsurf)
        logger.info(self,'Npang points %d' % self.npang)
        logger.info(self,'Ntrial %d' % self.ntrial)
        logger.info(self,'Rprimer %.6f' % self.rprimer)
        logger.debug(self, 'Rpru : %s' % self.rpru) 
        logger.info(self,'Epsiscp %.6f' % self.epsiscp)
        logger.info(self,'Epsroot %.6f' % self.epsroot)
        logger.info(self,'ODE solver %s' % self.backend)
        logger.info(self,'ODE tool %.6f' % self.epsilon)
        logger.info(self,'Max steps in ODE solver %d' % self.mstep)
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
            self.charges = f['molecule/charges'].value
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

        # 2) Setup basis info and grids
        if (self.ntrial%2 == 0): self.ntrial += 1
        geofac = numpy.power(((self.rmaxsurf-0.1)/self.rprimer),(1.0/(self.ntrial-1.0)))
        self.rpru = numpy.zeros((self.ntrial))
        for i in range(self.ntrial): 
            self.rpru[i] = self.rprimer*numpy.power(geofac,(i+1)-1)
        self.xnuc = numpy.asarray(self.coords[self.inuc])
        self.rsurf = numpy.zeros((self.npang,self.ntrial))
        self.nlimsurf = numpy.zeros((self.npang), dtype=numpy.int32)
        self.grids = grids.lebgrid(self.npang)
        
        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        # 3) Check rho nuclear atractors
        self.xyzrho = numpy.zeros((self.natm,3))
        for i in range(self.natm):
            self.xyzrho[i], gradmod = gradrho(self,self.coords[i],self.step)
            if (gradmod > 1e-4):
                if (self.charges[self.inuc] > 2.0):
                    logger.info(self,'Check rho position %.6f %.6f %.6f', *self.xyzrho[i])
                else:
                    raise RuntimeError('Failed finding nucleus:', *self.xyzrho[i]) 
            else:
                logger.info(self,'Check rho position %.6f %.6f %.6f', *self.xyzrho[i])
                logger.info(self,'Setting xyrho for atom to imput coords')
                self.xyzrho[i] = self.coords[i]

        backend = 1
        ct_ = numpy.asarray(self.grids[:,0], order='C')
        st_ = numpy.asarray(self.grids[:,1], order='C')
        cp_ = numpy.asarray(self.grids[:,2], order='C')
        sp_ = numpy.asarray(self.grids[:,3], order='C')
        angw_ = numpy.asarray(self.grids[:,4], order='C')
        # 4) Compute surface
        if (self.csurf):
            #pass
            feval = 'csurf_driver'
            drv = getattr(libcapi, feval)
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
                ctypes.c_int(self.npang),  
                ctypes.c_int((self.inuc)),  
                self.xyzrho.ctypes.data_as(ctypes.c_void_p), 
                ct_.ctypes.data_as(ctypes.c_void_p),
                st_.ctypes.data_as(ctypes.c_void_p),
                cp_.ctypes.data_as(ctypes.c_void_p),
                sp_.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(backend),
                ctypes.c_int(self.ntrial), 
                ctypes.c_double(self.epsiscp), 
                ctypes.c_double(self.epsroot), 
                ctypes.c_double(self.rmaxsurf), 
                ctypes.c_double(self.epsilon), 
                self.rpru.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_double(self.step), 
                ctypes.c_int(self.mstep),
                self.nlimsurf.ctypes.data_as(ctypes.c_void_p),
                self.rsurf.ctypes.data_as(ctypes.c_void_p))
        else:
            pass
        feval = 'surf_driver'
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
            ctypes.c_int(self.npang),  
            ctypes.c_int((self.inuc+1)),  
            self.xyzrho.ctypes.data_as(ctypes.c_void_p), 
            ctypes.c_char_p(self.chkfile),
            ct_.ctypes.data_as(ctypes.c_void_p),
            st_.ctypes.data_as(ctypes.c_void_p),
            cp_.ctypes.data_as(ctypes.c_void_p),
            sp_.ctypes.data_as(ctypes.c_void_p),
            angw_.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(backend),
            ctypes.c_int(self.ntrial), 
            ctypes.c_double(self.epsiscp), 
            ctypes.c_double(self.epsroot), 
            ctypes.c_double(self.rmaxsurf), 
            ctypes.c_double(self.epsilon), 
            ctypes.c_double(self.rprimer), 
            ctypes.c_double(self.step), 
            ctypes.c_int(self.mstep),
            self.nlimsurf.ctypes.data_as(ctypes.c_void_p),
            self.rsurf.ctypes.data_as(ctypes.c_void_p))

        self.rmin = 1000.0
        self.rmax = 0.0
        for i in range(self.npang):
            nsurf = int(self.nlimsurf[i])
            self.rmin = numpy.minimum(self.rmin,self.rsurf[i,0])
            self.rmax = numpy.maximum(self.rmax,self.rsurf[i,nsurf-1])
        logger.info(self,'Rmin for surface %.6f', self.rmin)
        logger.info(self,'Rmax for surface %.6f', self.rmax)

        # 5) Safe surface
        logger.info(self,'Finish with surface')
        logger.info(self,'Write HDF5 surface file')
        atom_dic = {'inuc':self.inuc,
                    'xnuc':self.xnuc,
                    'xyzrho':self.xyzrho,
                    'coords':self.grids,
                    'npang':self.npang,
                    'ntrial':self.ntrial,
                    'rmin':self.rmin,
                    'rmax':self.rmax,
                    'nlimsurf':self.nlimsurf,
                    'rsurf':self.rsurf}
        chkfile.save(self.chkfile, 'atom'+str(self.inuc), atom_dic)
        logger.info(self,'Surface of atom %d saved',self.inuc)
        logger.timer(self,'BaderSurf build', t0)
        logger.info(self,'')
    
        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    surf = BaderSurf(name)
    surf.epsilon = 1e-5
    surf.epsroot = 1e-5
    surf.verbose = 4
    surf.epsiscp = 0.220
    surf.mstep = 240
    surf.inuc = 0
    surf.npang = 14
    surf.csurf = True
    surf.kernel()
 
