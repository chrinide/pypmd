#!/usr/bin/env python

import os
import sys
import time
import h5py
import numpy
import ctypes
import signal
from pyscf import lib, dft
from pyscf.lib import logger

import grid

signal.signal(signal.SIGINT, signal.SIG_DFL)

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

_loaderpath = os.path.dirname(__file__)
libaim = numpy.ctypeslib.load_library('libaim.so', _loaderpath)

EPS = 1e-7

# TODO: screaning of points
def rho(self,x):
    x = numpy.reshape(x, (-1,3))
    ao = dft.numint.eval_ao(self.mol, x, deriv=0)
    ngrids, nao = ao.shape
    pos = self.mo_occ > self.occdrop
    cpos = numpy.einsum('ij,j->ij', self.mo_coeff[:,pos], numpy.sqrt(self.mo_occ[pos]))
    rho = numpy.zeros(ngrids)
    c0 = numpy.dot(ao, cpos)
    rho = numpy.einsum('pi,pi->p', c0, c0)
    return rho

def inbasin(self,r,j):
    isin = False
    rs1 = 0.0
    irange = self.nlimsurf[j]
    for k in range(irange):
        rs2 = self.rsurf[j,k]
        if (r >= rs1-EPS and r <= rs2+EPS):
            if (((k+1)%2) == 0):
                isin = False
            else:
                isin = True
            return isin
        rs1 = rs2
    return isin

def out_beta(self):
    logger.info(self,'* Go outside betasphere')
    xcoor = numpy.zeros(3)
    nrad = self.nrad
    npang = self.npang
    iqudr = self.iqudr
    mapr = self.mapr
    r0 = self.brad
    rfar = self.rmax
    rad = self.rad
    t0 = time.time()
    rmesh, rwei, dvol, dvoln = grid.rquad(nrad,r0,rfar,rad,iqudr,mapr)
    coordsang = self.agrids
    lmax = self.lmax
    NPROPS = lmax*(lmax+2) + 1
    t1 = time.time()
    ct_ = numpy.asarray(coordsang[:,0], order='C')
    st_ = numpy.asarray(coordsang[:,1], order='C')
    cp_ = numpy.asarray(coordsang[:,2], order='C')
    sp_ = numpy.asarray(coordsang[:,3], order='C')
    slm = numpy.zeros((npang,NPROPS))
    feval = 'eval_rsh'
    drv = getattr(libaim, feval)
    drv(ctypes.c_int(lmax), 
        ctypes.c_int(npang), 
        ct_.ctypes.data_as(ctypes.c_void_p),
        st_.ctypes.data_as(ctypes.c_void_p),
        cp_.ctypes.data_as(ctypes.c_void_p),
        sp_.ctypes.data_as(ctypes.c_void_p),
        slm.ctypes.data_as(ctypes.c_void_p)) 
    logger.debug(self,'Time finding outside RSH %.3f (sec)' % (time.time()-t1))
    rprops = numpy.zeros(NPROPS)
    for n in range(nrad):
        r = rmesh[n]
        coords = []
        weigths = []
        rslm = []
        for j in range(self.npang):
            inside = True
            inside = inbasin(self,r,j)
            if (inside == True):
                cost = coordsang[j,0]
                sintcosp = coordsang[j,1]*coordsang[j,2]
                sintsinp = coordsang[j,1]*coordsang[j,3]
                xcoor[0] = r*sintcosp
                xcoor[1] = r*sintsinp
                xcoor[2] = r*cost    
                p = self.xnuc + xcoor
                coords.append(p)
                rslm.append(slm[j,:])
                weigths.append(coordsang[j,4])
        coords = numpy.array(coords)
        weigths = numpy.array(weigths)
        rslm = numpy.array(rslm)
        val = numpy.einsum('i,ip->ip',rho(self,coords),rslm)
        props = numpy.einsum('ip,i->p', val, weigths)
        for l in range(lmax+1):
            ll = l*(l+1)
            for m in range(-l,l+1,1):
                lm = ll + m
                props[lm] *= r**l
        rprops += props*dvol[n]*rwei[n]
    logger.info(self,'*--> Qlm(0,0)  (s)      outside bsphere %f', rprops[0])    
    logger.info(self,'*--> Qlm(1,-1) (py)     outside bsphere %f', rprops[1])    
    logger.info(self,'*--> Qlm(1,0)  (pz)     outside bsphere %f', rprops[2])    
    logger.info(self,'*--> Qlm(1,1)  (px)     outside bsphere %f', rprops[3])    
    logger.info(self,'*--> Qlm(2,-2) (dxy)    outside bsphere %f', rprops[4])    
    logger.info(self,'*--> Qlm(2,-1) (dyz)    outside bsphere %f', rprops[5])    
    logger.info(self,'*--> Qlm(2,0)  (dz2)    outside bsphere %f', rprops[6])    
    logger.info(self,'*--> Qlm(2,1)  (dxz)    outside bsphere %f', rprops[7])    
    logger.info(self,'*--> Qlm(2,2)  (dx2-y2) outside bsphere %f', rprops[8])    
    logger.info(self,'Time out Bsphere %.3f (sec)' % (time.time()-t0))
    return rprops
    
def int_beta(self): 
    logger.info(self,'* Go with inside betasphere')
    xcoor = numpy.zeros(3)
    npang = self.bnpang
    coords = numpy.empty((npang,3))
    nrad = self.bnrad
    iqudr = self.biqudr
    mapr = self.bmapr
    r0 = 0
    rfar = self.brad
    rad = self.rad
    t0 = time.time()
    rmesh, rwei, dvol, dvoln = grid.rquad(nrad,r0,rfar,rad,iqudr,mapr)
    coordsang = grid.lebgrid(npang)
    lmax = self.blmax
    NPROPS = lmax*(lmax+2) + 1
    t1 = time.time()
    ct_ = numpy.asarray(coordsang[:,0], order='C')
    st_ = numpy.asarray(coordsang[:,1], order='C')
    cp_ = numpy.asarray(coordsang[:,2], order='C')
    sp_ = numpy.asarray(coordsang[:,3], order='C')
    slm = numpy.zeros((npang,NPROPS))
    feval = 'eval_rsh'
    drv = getattr(libaim, feval)
    drv(ctypes.c_int(lmax), 
        ctypes.c_int(npang), 
        ct_.ctypes.data_as(ctypes.c_void_p),
        st_.ctypes.data_as(ctypes.c_void_p),
        cp_.ctypes.data_as(ctypes.c_void_p),
        sp_.ctypes.data_as(ctypes.c_void_p),
        slm.ctypes.data_as(ctypes.c_void_p)) 
    logger.debug(self,'Time finding inside RSH %.3f (sec)' % (time.time()-t1))
    rprops = numpy.zeros(NPROPS)
    for n in range(nrad):
        r = rmesh[n]
        for j in range(npang): # j-loop can be changed to map
            cost = coordsang[j,0]
            sintcosp = coordsang[j,1]*coordsang[j,2]
            sintsinp = coordsang[j,1]*coordsang[j,3]
            xcoor[0] = r*sintcosp
            xcoor[1] = r*sintsinp
            xcoor[2] = r*cost    
            p = self.xnuc + xcoor
            coords[j] = p
        val = numpy.einsum('i,ip->ip',rho(self,coords),slm)
        props = numpy.einsum('ip,i->p', val, coordsang[:,4])
        for l in range(lmax+1):
            ll = l*(l+1)
            for m in range(-l,l+1,1):
                lm = ll + m
                props[lm] *= r**l
        rprops += props*dvol[n]*rwei[n]
    logger.info(self,'*--> Qlm(0,0)  (s)      inside bsphere %f', rprops[0])    
    logger.info(self,'*--> Qlm(1,-1) (py)     inside bsphere %f', rprops[1])    
    logger.info(self,'*--> Qlm(1,0)  (pz)     inside bsphere %f', rprops[2])    
    logger.info(self,'*--> Qlm(1,1)  (px)     inside bsphere %f', rprops[3])    
    logger.info(self,'*--> Qlm(2,-2) (dxy)    inside bsphere %f', rprops[4])    
    logger.info(self,'*--> Qlm(2,-1) (dyz)    inside bsphere %f', rprops[5])    
    logger.info(self,'*--> Qlm(2,0)  (dz2)    inside bsphere %f', rprops[6])    
    logger.info(self,'*--> Qlm(2,1)  (dxz)    inside bsphere %f', rprops[7])    
    logger.info(self,'*--> Qlm(2,2)  (dx2-y2) inside bsphere %f', rprops[8])    
    logger.info(self,'Time in Bsphere %.3f (sec)' % (time.time()-t0))
    return rprops

class Qlm(lib.StreamObject):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = lib.param.MAX_MEMORY
        self.chkfile = datafile
        self.surfile = datafile+'.h5'
        self.scratch = lib.param.TMPDIR 
        self.nthreads = lib.num_threads()
        self.inuc = 0
        self.nrad = 101
        self.iqudr = 'legendre'
        self.mapr = 'becke'
        self.betafac = 0.4
        self.bnrad = 101
        self.bnpang = 3074
        self.biqudr = 'legendre'
        self.bmapr = 'becke'
        self.non0tab = False
        self.corr = False
        self.occdrop = 1e-6
        self.lmax = 0
        self.blmax = 0
##################################################
# don't modify the following attributes, they are not input options
        self.rdm1 = None
        self.nocc = None
        self.mol = None
        self.mo_coeff = None
        self.mo_occ = None
        self.ntrial = None
        self.npang = None
        self.natm = None
        self.coords = None
        self.charges = None
        self.xnuc = None
        self.xyzrho = None
        self.agrids = None
        self.rsurf = None
        self.nlimsurf = None
        self.rmin = None
        self.rmax = None
        self.nelectron = None
        self.charge = None
        self.spin = None
        self.nprims = None
        self.nmo = None
        self.rad = None
        self.brad = None
        self._keys = set(self.__dict__.keys())

    def dump_input(self):

        if self.verbose < logger.INFO:
            return self

        logger.info(self,'')
        logger.info(self,'******** %s flags ********', self.__class__)
        logger.info(self,'* General Info')
        logger.info(self,'Date %s' % time.ctime())
        logger.info(self,'Python %s' % sys.version)
        logger.info(self,'Numpy %s' % numpy.__version__)
        logger.info(self,'Number of threads %d' % self.nthreads)
        logger.info(self,'Verbose level %d' % self.verbose)
        logger.info(self,'Scratch dir %s' % self.scratch)
        logger.info(self,'Input data file %s' % self.chkfile)
        logger.info(self,'Max_memory %d MB (current use %d MB)',
                 self.max_memory, lib.current_memory()[0])
        logger.info(self,'Correlated ? %s' % self.corr)

        logger.info(self,'* Molecular Info')
        logger.info(self,'Num atoms %d' % self.natm)
        logger.info(self,'Num electrons %d' % self.nelectron)
        logger.info(self,'Total charge %d' % self.charge)
        logger.info(self,'Spin %d ' % self.spin)
        logger.info(self,'Atom Coordinates (Bohr)')
        for i in range(self.natm):
            logger.info(self,'Nuclei %d with charge %d position : %.6f  %.6f  %.6f', i, 
                        self.charges[i], *self.coords[i])

        logger.info(self,'* Basis Info')
        logger.info(self,'Number of molecular orbitals %d' % self.nmo)
        logger.info(self,'Orbital EPS occ criterion %e' % self.occdrop)
        logger.info(self,'Number of occupied molecular orbitals %d' % self.nocc)
        logger.info(self,'Number of molecular primitives %d' % self.nprims)
        logger.debug(self,'Occs : %s' % self.mo_occ) 

        logger.info(self,'* Surface Info')
        logger.info(self,'Surface file %s' % self.surfile)
        logger.info(self,'Properties for nuc %d' % self.inuc)
        logger.info(self,'Nuclear coordinate %.6f  %.6f  %.6f', *self.xnuc)
        logger.info(self,'Rho nuclear coordinate %.6f  %.6f  %.6f', *self.xyzrho[self.inuc])
        logger.info(self,'Npang points %d' % self.npang)
        logger.info(self,'Ntrial %d' % self.ntrial)
        logger.info(self,'Rmin for surface %f', self.rmin)
        logger.info(self,'Rmax for surface %f', self.rmax)

        logger.info(self,'* Radial and angular grid Info')
        logger.info(self,'Npang points inside %d' % self.bnpang)
        logger.info(self,'Number of radial points outside %d', self.nrad)
        logger.info(self,'Number of radial points inside %d', self.bnrad)
        logger.info(self,'Radial outside quadrature %s', self.iqudr)
        logger.info(self,'Radial outside mapping %s', self.mapr)
        logger.info(self,'Radial inside quadrature %s', self.biqudr)
        logger.info(self,'Radial inside mapping %s', self.bmapr)
        logger.info(self,'Slater-Bragg radii %f', self.rad) 
        logger.info(self,'Beta-Sphere factor %f', self.betafac)
        logger.info(self,'Beta-Sphere radi %f', self.brad)

        logger.info(self,'* Real Spherical Harmonics expansion')
        logger.info(self,'Lmax inside %d' % self.blmax)
        logger.info(self,'Lmax outside %d', self.lmax)
        logger.info(self,'')

        return self

    def build(self):

        t0 = time.clock()
        lib.logger.TIMER_LEVEL = 3

        self.mol = lib.chkfile.load_mol(self.chkfile)
        self.nelectron = self.mol.nelectron 
        self.charge = self.mol.charge    
        self.spin = self.mol.spin      
        self.natm = self.mol.natm		
        self.coords = numpy.asarray([(numpy.asarray(atom[1])).tolist() for atom in self.mol._atom])
        self.charges = self.mol.atom_charges()
        self.mo_coeff = lib.chkfile.load(self.chkfile, 'scf/mo_coeff')
        self.mo_occ = lib.chkfile.load(self.chkfile, 'scf/mo_occ')
        nprims, nmo = self.mo_coeff.shape 
        self.nprims = nprims
        self.nmo = nmo
        if self.charges[self.inuc] == 1:
            self.rad = grid.BRAGG[self.charges[self.inuc]]
        else:
            self.rad = grid.BRAGG[self.charges[self.inuc]]*0.5

        if (self.corr):
            self.rdm1 = lib.chkfile.load(self.chkfile, 'rdm/rdm1') 
            natocc, natorb = numpy.linalg.eigh(self.rdm1)
            natorb = numpy.dot(self.mo_coeff, natorb)
            self.mo_coeff = natorb
            self.mo_occ = natocc
        nocc = self.mo_occ[abs(self.mo_occ)>self.occdrop]
        nocc = len(nocc)
        self.nocc = nocc

        idx = 'atom'+str(self.inuc)
        with h5py.File(self.surfile) as f:
            self.xnuc = f[idx+'/xnuc'].value
            self.xyzrho = f[idx+'/xyzrho'].value
            self.npang = f[idx+'/npang'].value
            self.ntrial = f[idx+'/ntrial'].value
            self.rmin = f[idx+'/rmin'].value
            self.rmax = f[idx+'/rmax'].value
            self.rsurf = f[idx+'/rsurf'].value
            self.nlimsurf = f[idx+'/nlimsurf'].value
            self.agrids = f[idx+'/coords'].value

        self.brad = self.rmin*self.betafac

        if self.verbose >= logger.WARN:
            self.check_sanity()
        if self.verbose > logger.NOTE:
            self.dump_input()

        if (self.iqudr == 'legendre'):
            self.iqudr = 1
        if (self.biqudr == 'legendre'):
            self.biqudr = 1

        if (self.mapr == 'becke'):
            self.mapr = 1
        elif (self.mapr == 'exp'):
            self.mapr = 2
        elif (self.mapr == 'none'):
            self.mapr = 0 
        if (self.bmapr == 'becke'):
            self.bmapr = 1
        elif (self.bmapr == 'exp'):
            self.bmapr = 2
        elif (self.bmapr == 'none'):
            self.bmapr = 0

        with lib.with_omp_threads(self.nthreads):
            brprops = int_beta(self)
            rprops = out_beta(self)

        logger.info(self,'Write info to HDF5 file')
        atom_dic = {'inprops':brprops,
                    'outprops':rprops,
                    'blmax':self.blmax,
                    'lmax':self.lmax,
                    'totprops':(brprops+rprops)}
        lib.chkfile.save(self.surfile, 'qlm'+str(self.inuc), atom_dic)

        logger.info(self,'*-> Total Qlm(0,0)  (s)      %f' % (rprops[0]+brprops[0]))   
        logger.info(self,'*-> Total Qlm(1,-1) (py)     %f' % (rprops[1]+brprops[1]))   
        logger.info(self,'*-> Total Qlm(1,0)  (pz)     %f' % (rprops[2]+brprops[2]))   
        logger.info(self,'*-> Total Qlm(1,1)  (px)     %f' % (rprops[3]+brprops[3]))   
        logger.info(self,'*-> Total Qlm(2,-2) (dxy)    %f' % (rprops[4]+brprops[4]))   
        logger.info(self,'*-> Total Qlm(2,-1) (dyz)    %f' % (rprops[5]+brprops[5]))   
        logger.info(self,'*-> Total Qlm(2,0)  (dz2)    %f' % (rprops[6]+brprops[6]))   
        logger.info(self,'*-> Total Qlm(2,1)  (xz)     %f' % (rprops[7]+brprops[7]))   
        logger.info(self,'*-> Total Qlm(2,2)  (dx2-y2) %f' % (rprops[8]+brprops[8]))   
        logger.info(self,'')

        logger.info(self,'Qlm of atom %d done',self.inuc)
        logger.timer(self,'Qlm build', t0)

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.chk'
    natm = 3
    bas = Qlm(name)
    bas.verbose = 4
    bas.nrad = 221
    bas.iqudr = 'legendre'
    bas.mapr = 'becke'
    bas.bnrad = 121
    bas.bnpang = 3074
    bas.biqudr = 'legendre'
    bas.bmapr = 'becke'
    bas.betafac = 0.4
    bas.non0tab = False
    bas.lmax = 10
    bas.blmax = 10
    for i in range(natm):
        bas.inuc = i
        bas.kernel()

