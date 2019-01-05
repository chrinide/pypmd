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
import numint

signal.signal(signal.SIGINT, signal.SIG_DFL)

libfapi = misc.load_library('libfapi')

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

EPS = 1e-7

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

# TODO: better iqudr and mapr selection
def out_beta(self):
    logger.info(self,'* Go outside betasphere')
    xcoor = numpy.zeros(3)
    nrad = self.nrad
    iqudr = self.iqudr
    mapr = self.mapr
    r0 = self.brad
    rfar = self.rmax
    rad = self.rad
    t0 = time.clock()
    rmesh, rwei, dvol, dvoln = grids.rquad(nrad,r0,rfar,rad,iqudr,mapr)
    coordsang = self.agrids
    rprops = 0.0
    for n in range(nrad):
        r = rmesh[n]
        coords = []
        weigths = []
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
                p = self.xyzrho + xcoor
                coords.append(p)
                weigths.append(coordsang[j,4])
        coords = numpy.array(coords)
        weigths = numpy.array(weigths)
        val = numint.eval_rho(self,coords)
        props = numpy.einsum('i,i->', val, weigths)
        rprops += props*dvol[n]*rwei[n]
    logger.info(self,'*--> Electron density outside bsphere %8.5f ', rprops)    
    logger.timer(self,'Out Bsphere build', t0)
    return rprops
    
# TODO: better iqudr and mapr selection
def int_beta(self): 
    logger.info(self,'* Go with inside betasphere')
    xcoor = numpy.zeros(3)
    coords = numpy.empty((self.bnpang,3))
    nrad = self.bnrad
    iqudr = self.biqudr
    mapr = self.bmapr
    r0 = 0
    rfar = self.brad
    rad = self.rad
    t0 = time.clock()
    rmesh, rwei, dvol, dvoln = grids.rquad(nrad,r0,rfar,rad,iqudr,mapr)
    coordsang = grids.lebgrid(self.bnpang)
    rprops = 0.0
    for n in range(nrad):
        r = rmesh[n]
        for j in range(self.bnpang): # j-loop can be changed to map
            cost = coordsang[j,0]
            sintcosp = coordsang[j,1]*coordsang[j,2]
            sintsinp = coordsang[j,1]*coordsang[j,3]
            xcoor[0] = r*sintcosp
            xcoor[1] = r*sintsinp
            xcoor[2] = r*cost    
            p = self.xnuc + xcoor
            coords[j] = p
        val = numint.eval_rho(self,coords)
        props = numpy.einsum('i,i->', val, coordsang[:,4])
        rprops += props*dvol[n]*rwei[n]
    logger.info(self,'*--> Electron density inside bsphere %8.5f ', rprops)    
    logger.timer(self,'Bsphere build', t0)
    return rprops

class Basin(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
        self.small_rho_cutoff = 1e-7
        self.nthreads = misc.num_threads()
        self.gpu = False
        self.inuc = 0
        self.nrad = 101
        self.iqudr = 'legendre'
        self.mapr = 'basin'
        self.betafac = 0.4
        self.bnrad = 101
        self.biqudr = 'legendre'
        self.bmapr = 'basin'
        self.bnpang = 3074
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
        self.ntrial = None
        self.npang = None
        self.xnuc = None
        self.xyzrho = None
        self.agrids = None
        self.rsurf = None
        self.nlimsurf = None
        self.rmin = None
        self.rmax = None
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

        logger.info(self,'* Surface Info')
        logger.info(self,'Comuting properties for atom %d' % self.inuc)
        logger.info(self,'Nuclei with position : %.6f  %.6f  %.6f', *self.xnuc)
        logger.info(self,'Number of angular points %d' % self.npang)
        logger.info(self,'Rmin for surface %.6f', self.rmin) 
        logger.info(self,'Rmax for surface %.6f', self.rmax) 

        logger.info(self,'* Integration grids Info')
        logger.info(self,'Number of radial points outside %d', self.nrad)
        logger.info(self,'Number of radial points inside %d', self.bnrad)
        logger.info(self,'Radial outside quadrature %s', self.iqudr)
        logger.info(self,'Radial outside mapping %s', self.mapr)
        logger.info(self,'Radial inside quadrature %s', self.biqudr)
        logger.info(self,'Radial inside mapping %s', self.bmapr)
        logger.info(self,'Slater-Bragg radii %.6f', self.rad) 
        logger.info(self,'Beta-Sphere factor %.6f', self.betafac)
        logger.info(self,'Beta-Sphere radi %.6f', self.brad)
        logger.info(self,'Number of angular points in beta %d' % self.bnpang)
        logger.info(self,'')

        return self

    def build(self):

        t0 = time.clock()
        logger.TIMER_LEVEL = 3
    
        # 1) Build info
        idx = 'atom'+str(self.inuc)
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
            self.inuc = f[idx+'/inuc'].value
            self.xnuc = f[idx+'/xnuc'].value
            self.xyzrho = f[idx+'/xyzrho'].value
            self.agrids = f[idx+'/coords'].value
            self.npang = f[idx+'/npang'].value
            self.ntrial = f[idx+'/ntrial'].value
            self.rmin = f[idx+'/rmin'].value
            self.rmax = f[idx+'/rmax'].value
            self.nlimsurf = f[idx+'/nlimsurf'].value
            self.rsurf = f[idx+'/rsurf'].value

        if self.charges[self.inuc] == 1:
            self.rad = grids.BRAGG[self.charges[self.inuc]]
        else:
            self.rad = grids.BRAGG[self.charges[self.inuc]]*0.5
        self.brad = self.rmin*self.betafac

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

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        # 2) Compute integrals
        brprops = int_beta(self)
        rprops = out_beta(self)

        #logger.info(self,'Write info to HDF5 file')
        #atom_dic = {'inprops':brprops,
        #            'outprops':rprops,
        #            'totprops':(brprops+rprops)}
        #chkfile.save(self.surfile, 'atom_props'+str(self.inuc), atom_dic)
        logger.info(self,'*-> Total density density %8.5f ', (rprops+brprops))    
        logger.info(self,'')
        logger.info(self,'Basim properties of atom %d done',self.inuc)
        logger.timer(self,'Basin build', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    bas = Basin(name)
    bas.verbose = 4
    bas.nrad = 101
    bas.iqudr = 'legendre'
    bas.mapr = 'becke'
    bas.bnrad = 101
    bas.bnpang = 5810
    bas.biqudr = 'legendre'
    bas.bmapr = 'becke'

    bas.inuc = 0
    bas.kernel()

    bas.inuc = 1
    bas.kernel()
 
    #bas.inuc = 2
    #bas.kernel()

