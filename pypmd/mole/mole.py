#!/usr/bin/env python

import os
import sys
import time
import numpy
import ctypes

from pypmd import lib
from pypmd.lib import logger
from pypmd import io

libwfn = lib.load_library('libsurf')

MGRP = 500
NGTOH = 21
MAXTYPE = 56

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

class Mole(lib.StreamObject):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = lib.param.MAX_MEMORY
        self.wfnfile = datafile
        self.scratch = lib.param.TMPDIR 
        self.epscuttz = 1e-14
        self.rmaxatom = 20.0
        self.epsortho = 1e-5
        self.check = True
##################################################
# don't modify the following attributes, they are not input options
        self.nmo = None
        self.mo_count = None
        self.nprim = None
        self.atnam = None
        self.icen = None
        self.ityp = None
        self.oexp = None
        self.mo_coeff = None
        self.mo_occ = None
        self.mo_energy = None
        self.natm = None
        self.coords = None
        self.charges = None
        self.nelectrons = None
        self.charge = None
        self.spin = None
        self.cart = None
        self._built = None
        self._keys = set(self.__dict__.keys())

    def dump_input(self):

        if self.verbose < logger.INFO:
            return self

        logger.info(self,'')
        logger.info(self,'******** %s flags ********', self.__class__)
        logger.info(self,'Verbose level %d' % self.verbose)
        logger.info(self,'Scratch dir %s' % self.scratch)
        logger.info(self,'Input data file %s' % self.wfnfile)
        logger.info(self,'max_memory %d MB ', self.max_memory)
        logger.info(self,'Num atoms %d' % self.natm)
        logger.info(self,'Num electrons %d' % self.nelectrons)
        logger.info(self,'Total charge %d' % self.charge)
        logger.info(self,'Atom Coordinates (Bohr)')
        for i in range(self.natm):
            logger.info(self,'Nuclei %d %s with charge %d position : %.6f  %.6f  %.6f', \
            i, self.atnam[i], self.charges[i], *self.coords[i])
        logger.info(self,'Total number of MO %d' % self.nmo)
        logger.info(self,'Total number of primitivies %d' % self.nprim)
        logger.info(self,'Screaning for primitivies %.6f' % self.epscuttz)
        logger.info(self,'Max distance for primitivie %.6f' % self.rmaxatom)
        logger.info(self,'Check wave function %s' % self.check)
        if (self.check):
            logger.info(self,'EPS for orthogonalization %.6f' % self.epsortho)
        return self

    def build(self):

        t0 = time.clock()
        lib.logger.TIMER_LEVEL = 5

        title, symbols, numbers, coords, icen, ityp, oexp, \
        mo_count, mo_occ, mo_energy, mo_coeff = io.load_wfn(self.wfnfile)
        nprims, nmo = mo_coeff.shape
        natm = coords.shape[0]

        self.atnam = symbols
        self.charges = numbers
        self.coords = numpy.asarray(coords, order='C')
        self.icen = icen
        self.ityp = ityp
        self.oexp = oexp
        self.mo_count = mo_count
        self.mo_coeff = numpy.asarray(mo_coeff, order='C')
        self.mo_occ = mo_occ
        self.mo_energy = mo_energy
        self.nmo = nmo
        self.nprim = nprims
        self.natm = natm
        self.nelectrons = mo_occ.sum()
        self.charge = numbers.sum() - self.nelectrons
        self.cart = True
        self._built = True

        if self.verbose >= logger.WARN:
            self.check_sanity()
        if self.verbose > lib.logger.NOTE:
            self.dump_input()

        newnprim = numpy.zeros(1, dtype=numpy.int32)
        numshells = numpy.zeros(1, dtype=numpy.int32)
        maxgrp = numpy.zeros(1, dtype=numpy.int32)
        npc = numpy.zeros(self.natm, dtype=numpy.int32)
        ngroup = numpy.zeros(self.natm, dtype=numpy.int32)
        icenat = numpy.zeros((self.nprim,self.natm), dtype=numpy.int32)
        nzexp = numpy.zeros((MGRP,self.natm), dtype=numpy.int32) 
        nuexp = numpy.zeros((NGTOH,MGRP,self.natm), dtype=numpy.int32)  
        rcutte = numpy.zeros((MGRP,self.natm), dtype=numpy.float64)  
    
        libwfn.wfn_driver(ctypes.c_double(self.epscuttz),  
                          ctypes.c_double(self.rmaxatom),  
                          ctypes.c_int(self.nmo),  
                          ctypes.c_int(self.nprim),  
                          self.icen.ctypes.data_as(ctypes.c_void_p), 
                          self.ityp.ctypes.data_as(ctypes.c_void_p), 
                          self.oexp.ctypes.data_as(ctypes.c_void_p), 
                          self.mo_coeff.ctypes.data_as(ctypes.c_void_p), 
                          self.mo_occ.ctypes.data_as(ctypes.c_void_p), 
                          ctypes.c_int(self.natm),  
                          self.coords.ctypes.data_as(ctypes.c_void_p), 
                          self.charges.ctypes.data_as(ctypes.c_void_p),
                          newnprim.ctypes.data_as(ctypes.c_void_p),
                          maxgrp.ctypes.data_as(ctypes.c_void_p),
                          numshells.ctypes.data_as(ctypes.c_void_p),
                          npc.ctypes.data_as(ctypes.c_void_p),
                          ngroup.ctypes.data_as(ctypes.c_void_p),
                          icenat.ctypes.data_as(ctypes.c_void_p),
                          nzexp.ctypes.data_as(ctypes.c_void_p),
                          nuexp.ctypes.data_as(ctypes.c_void_p),
                          rcutte.ctypes.data_as(ctypes.c_void_p))

        logger.timer(self,'Mol and wfn info build', t0)
        #print newnprim, maxgrp, numshells
        #print "npc", npc
        #print "ngroup", ngroup
        #for j in range(self.natm):
        #    for i in range(npc[j]):
        #        print icenat[i,j]
        #for j in range(self.natm):
        #    for i in range(ngroup[j]):
        #        print nzexp[i,j]
        #for j in range(self.natm):
        #    for i in range(ngroup[j]):
        #        for k in range(nzexp[i,j]):
        #            print nuexp[k,i,j]
        for j in range(self.natm):
            for i in range(ngroup[j]):
                print rcutte[i,j]

        return self

    kernel = build

if __name__ == '__main__':
    name = 'lif.wfn'
    mol = Mole(name)
    mol.verbose = 4
    mol.kernel()
 
    #integer(kind=ip), parameter :: ncentw = 100
    #character(len=1) :: lb(0:5), lbl
    #logical :: wrout
    #data (lb(i),i=0,5) /'S','P','D','F','G','H'/
    #if (ncent.gt.ncentw) then
    #  wrout = .false.
    #else
    #  wrout = .true.
    #end if
    #write (uout,'(1x,a,1x,i0)') string('# Actual number of primitives :'), nprims
    #write (uout,'(1x,a,1x,i0)') string('# Total number of shells :'), numshells
    #if (wrout) then
    #  do ic = 1,ncent
    #    write (uout,210) ic
    #    write (uout,300) nshell(ic),ic
    #    write (uout,301) (ishell(ic,j),atcenter(ic,j),j=1,nshell(ic))
    #    do m = 1,ngroup(ic)
    #      i = nuexp(ic,m,1)
    #      lsum = nlm(ityp(i),1)+nlm(ityp(i),2)+nlm(ityp(i),3)
    #      lbl = lb(lsum)
    #      write (uout,613) lbl,oexp(i),sqrt(rcutte(ic,m)),(nuexp(ic,m,k),k=1,nzexp(ic,m))
    #    end do
    #  end do
    #end if

#300 format (1x,'# ',i0,' shells contribute to the basin of center ',i0, &
#     /,' # [ shell(atom) means shell number "shell" of atom "atom" ]')
#301 format (1000(1x,8(I6,'(',I4,')'),/))
#210 format (1x,'# CENTER ',i0)
#613 format (1x,'# ',a,' Shell   Exp = ',e16.8,4x,'Cutoff = ',f13.6, &
#            4x,'Primitives : ',21(1x,i0))

