#!/usr/bin/env python

import sys
import h5py
import time
import numpy
import ctypes

import misc
import data
import param
import logger


# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str
    

libcapi = misc.load_library('libcapi.so')
libfapi = misc.load_library('libfapi.so')
        

def vv10(self,rho,gnorm):        
    gnorm2 = numpy.zeros(self.npoints)
    for i in range(self.npoints):
        gnorm2[i] = gnorm[i]*gnorm[i]
    coef_beta = 1.0/32.0 * (3.0/(self.coef_B**2.0))**(3.0/4.0)
    kappa_pref = self.coef_B * (1.5*numpy.pi)/((9.0*numpy.pi)**(1.0/6.0))
    const = 4.0/3.0 * numpy.pi
    idx = numpy.nonzero(rho)
    rho = rho[idx]
    gnorm2 = gnorm2[idx]
    ngrids = rho.shape[0]
    coords = self.points[idx]
    weights = self.weights[idx]
    vv10_e = 0.0
    for idx1 in range(ngrids):
        point1 = coords[idx1,:]
        rho1 = rho[idx1]
        weigth1 = weights[idx1]
        gamma1 = gnorm2[idx1]
        Wp1 = const*rho1
        Wg1 = self.coef_C * ((gamma1/(rho1*rho1))**2.0)
        W01 = numpy.sqrt(Wg1 + Wp1)
        kappa1 = rho1**(1.0/6.0)*kappa_pref
        R =  (point1[0]-coords[:,0])**2
        R += (point1[1]-coords[:,1])**2
        R += (point1[2]-coords[:,2])**2
        Wp2 = const*rho
        Wg2 = self.coef_C * ((gnorm2/(rho*rho))**2.0)
        W02 = numpy.sqrt(Wg2 + Wp2)
        kappa2 = rho**(1.0/6.0)*kappa_pref
        g = W01*R + kappa1
        gp = W02*R + kappa2
        kernel12 = -1.5*weights*rho/(g*gp*(g+gp))
        kernel = kernel12.sum()
        vv10_e += weigth1*rho1*(coef_beta + 0.5*kernel)
    return vv10_e


def c_vv10(self,rho,gnorm):        
    gnorm2 = numpy.zeros(self.npoints)
    for i in range(self.npoints):
        gnorm2[i] = gnorm[i]*gnorm[i]
    coef_beta = 1.0/32.0 * (3.0/(self.coef_B**2.0))**(3.0/4.0)
    kappa_pref = self.coef_B * (1.5*numpy.pi)/((9.0*numpy.pi)**(1.0/6.0))
    const = 4.0/3.0 * numpy.pi
    idx = numpy.nonzero(rho)
    rho = rho[idx]
    gnorm2 = gnorm2[idx]
    ngrids = rho.shape[0]
    coords = self.points[idx]
    weights = self.weights[idx]
    libcapi.vv10.restype = ctypes.c_double
    ev = libcapi.vv10(ctypes.c_int(ngrids),
             ctypes.c_double(self.coef_C),
             ctypes.c_double(self.coef_B),
             coords.ctypes.data_as(ctypes.c_void_p),
             rho.ctypes.data_as(ctypes.c_void_p),
             weights.ctypes.data_as(ctypes.c_void_p),
             gnorm2.ctypes.data_as(ctypes.c_void_p))
    return ev

    
class Disp(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
        self.coef_C = 0.0093
        self.coef_B = 5.9
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
        self.npoints = None
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
        logger.info(self,'Use GPU %s' % self.gpu)
        if (self.gpu):
            logger.info(self,'* GPU Info')

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

        logger.info(self,'* Grid and VV10 Info')
        logger.info(self,'Tot grids %d', self.npoints)
        logger.info(self,'Coef C %f', self.coef_C)
        logger.info(self,'Coef B %f', self.coef_B)
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
            self.npoints = f['becke/npoints'].value
            self.weights = f['becke/weights'].value
            self.points = f['becke/points'].value

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        # 2) Compute VV10 dispersion energy
        rho = numpy.zeros(self.npoints)
        gnorm = numpy.zeros(self.npoints)
        if (self.gpu):
            pass
        else:
            t1 = time.clock()
            libcapi.driver_vv10(ctypes.c_int(self.nmo),  
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
                           ctypes.c_int(self.npoints),  
                           self.points.ctypes.data_as(ctypes.c_void_p), 
                           rho.ctypes.data_as(ctypes.c_void_p), 
                           gnorm.ctypes.data_as(ctypes.c_void_p)) 
            logger.info(self,'Integral of Rho %.12f' % numpy.einsum('i,i->', rho, self.weights))
            logger.timer(self,'Rho and norm time', t1)
            t1 = time.clock()
            vv10_e = vv10(self,rho,gnorm)         
            logger.info(self, 'Python VV10 energy %.12f' % vv10_e)
            logger.timer(self,'Python time', t1)
            t1 = time.clock()
            vv10_e = c_vv10(self,rho,gnorm)         
            logger.info(self, 'C VV10 energy %.12f' % vv10_e)
            logger.timer(self,'C time', t1)

        logger.info(self,'Finish with vv10')
        logger.timer(self,'VV10 time', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    wfn = Disp(name)
    wfn.verbose = 4
    wfn.kernel()
