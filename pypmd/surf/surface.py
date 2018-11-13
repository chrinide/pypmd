#!/usr/bin/env python

import os
import sys
import time
import numpy
import ctypes

from pypmd import lib
from pypmd.lib import logger
from pypmd import grid

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

class BaderSurf(lib.StreamObject):

    def __init__(self, mol):
        self.mol = mol
        self.verbose = mol.verbose
        self.max_memory = mol.max_memory
        self.stdout = mol.stdout
        self.scratch = lib.param.TMPDIR 
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
        self.mstep = 100
        self.wfnfile = mol.wfnfile
        self.surfile = self.wfnfile+'.h5'
##################################################
# don't modify the following attributes, they are not input options
        self.mo_coeff = mol.mo_coeff
        self.mo_occ = mol.mo_occ
        self.natm = mol.natm
        self.nmo = mol.nmo
        self.coords = mol.coords
        self.charges = mol.charges
        self.nprim = mol.nprim
        self.cart = mol.cart
        self.xnuc = None
        self.xyzrho = None
        self.rpru = None
        self.grids = None
        self.rsurf = None
        self.nlimsurf = None
        self.rmin = None
        self.rmax = None
        self._keys = set(self.__dict__.keys())

    def dump_input(self):

        if self.verbose < logger.INFO:
            return self

        logger.info(self,'')
        logger.info(self,'******** %s flags ********', self.__class__)
        logger.info(self,'Verbose level %d' % self.verbose)
        logger.info(self,'Scratch dir %s' % self.scratch)
        logger.info(self,'Input data file %s' % self.chkfile)
        logger.info(self,'max_memory %d MB' % self.max_memory)
        logger.info(self,'Atom Coordinates (Bohr)')
        logger.info(self,'Surface file %s' % self.surfile)
        logger.info(self,'Surface for nuc %d' % self.inuc)
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

        return self

    def build(self):

        t0 = time.clock()
        lib.logger.TIMER_LEVEL = 5

        if (self.ntrial%2 == 0): self.ntrial += 1
        geofac = numpy.power(((self.rmaxsurf-0.1)/self.rprimer),(1.0/(self.ntrial-1.0)))
        self.rpru = numpy.zeros((self.ntrial))
        for i in range(self.ntrial): 
            self.rpru[i] = self.rprimer*numpy.power(geofac,(i+1)-1)
        self.xnuc = numpy.asarray(self.coords[self.inuc])
        self.rsurf = numpy.zeros((self.npang,self.ntrial))
        self.nlimsurf = numpy.zeros((self.npang), dtype=numpy.int32)
        #self.grids = grid.lebgrid(self.npang)

        if self.verbose >= logger.WARN:
            self.check_sanity()
        if self.verbose > lib.logger.NOTE:
            self.dump_input()

        self.xyzrho = numpy.zeros((3))
        #self.xyzrho, gradmod = cp.gradrho(self,self.xnuc,self.step)
        #if (gradmod > 1e-4):
        #    if (self.charges[self.inuc] > 2.0):
        #        logger.info(self,'Check rho position %.6f %.6f %.6f', *self.xyzrho)
        #    else:
        #        raise RuntimeError('Failed finding nucleus:', *self.xyzrho) 
        #else:
        #    logger.info(self,'Check rho position %.6f %.6f %.6f', *self.xyzrho)
        #
        # TODO: check why data is not persistent before C call
        #atom_dic = {'inuc':self.inuc,
        #            'xnuc':self.xnuc,
        #            'xyzrho':self.xyzrho,
        #            'coords':self.grids,
        #            'npang':self.npang,
        #            'ntrial':self.ntrial}
        #lib.chkfile.save(self.surfile, 'atom'+str(self.inuc), atom_dic)
        #
        #feval = 'surf_driver'
        #drv = getattr(libaim, feval)
        #backend = 1
        #ct_ = numpy.asarray(self.grids[:,0], order='C')
        #st_ = numpy.asarray(self.grids[:,1], order='C')
        #cp_ = numpy.asarray(self.grids[:,2], order='C')
        #sp_ = numpy.asarray(self.grids[:,3], order='C')
        #drv(ctypes.c_int(self.inuc), 
        #    ctypes.c_int(self.npang), 
        #    ct_.ctypes.data_as(ctypes.c_void_p),
        #    st_.ctypes.data_as(ctypes.c_void_p),
        #    cp_.ctypes.data_as(ctypes.c_void_p),
        #    sp_.ctypes.data_as(ctypes.c_void_p),
        #    ctypes.c_int(self.ntrial), self.rpru.ctypes.data_as(ctypes.c_void_p), 
        #    ctypes.c_double(self.epsiscp), ctypes.c_double(self.epsroot), 
        #    ctypes.c_double(self.rmaxsurf), ctypes.c_int(backend),
        #    ctypes.c_double(self.epsilon), ctypes.c_double(self.step), 
        #    ctypes.c_int(self.mstep),
        #    ctypes.c_int(self.cart),
        #    self.coords.ctypes.data_as(ctypes.c_void_p),
        #    self.xyzrho.ctypes.data_as(ctypes.c_void_p),
        #    self.atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(self.natm),
        #    self.bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(self.nbas),
        #    self.env.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(self.nprim),
        #    self.ao_loc.ctypes.data_as(ctypes.c_void_p),
        #    self.mo_coeff.ctypes.data_as(ctypes.c_void_p),
        #    self.mo_occ.ctypes.data_as(ctypes.c_void_p),
        #    self.nlimsurf.ctypes.data_as(ctypes.c_void_p),
        #    self.rsurf.ctypes.data_as(ctypes.c_void_p))
        #
        #self.rmin = 1000.0
        #self.rmax = 0.0
        #for i in range(self.npang):
        #    nsurf = int(self.nlimsurf[i])
        #    self.rmin = numpy.minimum(self.rmin,self.rsurf[i,0])
        #    self.rmax = numpy.maximum(self.rmax,self.rsurf[i,nsurf-1])
        #logger.info(self,'Rmin for surface %.6f', self.rmin)
        #logger.info(self,'Rmax for surface %.6f', self.rmax)
        #
        #logger.info(self,'Write HDF5 surface file')
        #lib.chkfile.save(self.surfile, 'atom'+str(self.inuc)+'/intersecs', self.nlimsurf)
        #lib.chkfile.save(self.surfile, 'atom'+str(self.inuc)+'/surface', self.rsurf)
        #lib.chkfile.save(self.surfile, 'atom'+str(self.inuc)+'/rmin', self.rmin)
        #lib.chkfile.save(self.surfile, 'atom'+str(self.inuc)+'/rmax', self.rmax)
        #logger.info(self,'Surface of atom %d saved',self.inuc)
        #logger.timer(self,'BaderSurf build', t0)
        #
        return self

    kernel = build

if __name__ == '__main__':

    from pypmd.mole import Mole    
    name = 'lif.wfn'
    mol = Mole(name)
    mol.verbose = 4
    mol.kernel()

    surf = BaderSurf(mol)
    surf.epsilon = 1e-4
    surf.epsiscp = 0.180
    surf.mstep = 100
    surf.inuc = 0
    surf.npang = 5810
    surf.kernel()
 
