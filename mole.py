#!/usr/bin/env python

import os
import sys
import time
import numpy
import ctypes

import wfn
import data
import param
import misc
import chkfile
import logger
import evaluator
import tools

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str


def purge_prims(self):    
    self.oldnprims = self.nprims
    keep = numpy.ones(self.nprims)
    for i in range(self.nprims):
        if (keep[i] == 0): continue
        for j in range(i):
            if (abs(self.oexp[j] - self.oexp[i]) >= 1e-10): continue
            if (self.ityp[i] == self.ityp[j] & self.icen[i] == self.icen[j]):
                #mo_coeff[i,:] += mo_coeff[j,:]
                keep[j] = 0
    #print "Keep", keep
    return self


def get_shells(self):
    icenat = numpy.zeros((self.nprims, self.natm), dtype=numpy.int32)
    for i in range(self.nprims):
        itip = self.ityp[i] - 1
        it0 = data.nlm[itip, 0]
        it1 = data.nlm[itip, 1]
        it2 = data.nlm[itip, 2]
        self.lang[i] = it0 + it1 + it2
        ic = self.icen[i] - 1
        self.npc[ic] += 1
        inda = self.npc[ic] - 1
        icenat[inda, ic] = i + 1
    self.lmax = max(self.lang)
    self.ngroup = numpy.zeros(self.natm, dtype=numpy.int32)
    self.nzexp = numpy.zeros((param.MGRP, self.natm), dtype=numpy.int32)
    self.nuexp = numpy.zeros((param.MGRP, param.NGTOH, self.natm), dtype=numpy.int32)
    for ic in range(self.natm):
        for j in range(self.npc[ic]):
            cond = False
            k = icenat[j, ic]
            itk = self.ityp[k - 1] - 1
            it0 = data.nlm[itk, 0]
            it1 = data.nlm[itk, 1]
            it2 = data.nlm[itk, 2]
            isuk = it0 + it1 + it2
            if (j == 0):
                self.ngroup[ic] = 1
                self.nzexp[0, ic] = 1
                self.nuexp[0, 0, ic] = k
            else:
                for m in range(self.ngroup[ic]):
                    inda = self.nuexp[m, 0, ic]
                    itd = self.ityp[inda - 1] - 1
                    it0 = data.nlm[itd, 0]
                    it1 = data.nlm[itd, 1]
                    it2 = data.nlm[itd, 2]
                    isud = it0 + it1 + it2
                    if (abs(self.oexp[k - 1] - self.oexp[inda - 1]) <= 1e-8):
                        if (itk == 1 & itd == 1):
                            #raise RuntimeError('''Two S primitives with equal exponents''')
                            pass
                        else:
                            if (isuk == isud):
                                self.nzexp[m, ic] += 1
                                self.nuexp[m, self.nzexp[m, ic] - 1, ic] = k
                                cond = True
                                break
                if (cond != True):
                    self.ngroup[ic] += 1
                    if (self.ngroup[ic] >= param.MGRP):
                        raise RuntimeError('''Increase MGRP''')
                    self.nzexp[self.ngroup[ic] - 1, ic] = 1
                    self.nuexp[self.ngroup[ic] - 1, 0, ic] = k

    self.nshells = self.ngroup.sum() 
    return self


def get_shells_eps(self): 
    self.rcutte = numpy.zeros((param.MGRP,self.natm))
    for ic in range(self.natm):
        for m in range(self.ngroup[ic]):
          i = self.nuexp[m,0,ic]
          itip = self.ityp[i - 1] - 1
          it1 = data.nlm[itip, 0]
          it2 = data.nlm[itip, 1]
          it3 = data.nlm[itip, 2]
          lsum = it1 + it2 + it3
          zz = self.oexp[i-1]
          x1 = 0.1
          while (x1**lsum*numpy.exp(-zz*x1*x1)>=numpy.abs(self.cuttz)):
            x1 += 0.1
          self.rcutte[m,ic] = x1*x1
    return self


class Mole(object):

    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.wfnfile = datafile
        self.scratch = param.TMPDIR 
        self.cuttz = 1e-12
        self.chkfile = datafile+'.h5'
        self.rdmfile = datafile+'.rdm1'
        self.nthreads = misc.num_threads()
##################################################
# don't modify the following attributes, they are not input options
        self.natm = None
        self.symbols = None
        self.coords = None
        self.charges = None
        self.nelectrons = None
        self.charge = None
        self.spin = None
        self.nmo = None
        self.nprims = None
        self.oldnprims = None
        self.icen = None
        self.ityp = None
        self.oexp = None
        self.lang = None
        self.lmax = None
        self.npc = None
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
        logger.info(self,'Date %s' % time.ctime())
        logger.info(self,'Python %s' % sys.version)
        logger.info(self,'Numpy %s' % numpy.__version__)
        logger.info(self,'Number of threads %d' % self.nthreads)
        logger.info(self,'Verbose level %d' % self.verbose)
        logger.info(self,'Scratch dir %s' % self.scratch)
        logger.info(self,'Input chkfile data file %s' % self.chkfile)
        logger.info(self,'Max_memory %d MB' % self.max_memory)
        logger.info(self,'Input wfn data file %s' % self.wfnfile)
        logger.info(self,'Output h5 data file %s' % self.chkfile)

        logger.info(self,'* Molecular Info')
        logger.info(self,'Num atoms %d' % self.natm)
        logger.info(self,'Num electrons %d' % self.nelectrons)
        #logger.info(self,'Total charge %d' % self.charge)
        #logger.info(self,'Spin %d ' % self.spin)
        logger.info(self,'Atom Coordinates (Bohr)')
        for i in range(self.natm):
            logger.info(self,'Nuclei %d %s with charge %d position : %.6f  %.6f  %.6f', i, 
                        self.symbols[i], self.charges[i], *self.coords[i])

        logger.info(self,'* Basis Info')
        logger.info(self,'Cuttoff for primitives %g' % self.cuttz)
        logger.info(self,'Number of Orbitals %d' % self.nmo)
        logger.info(self,'Old Number of primitives %d' % self.oldnprims)
        logger.info(self,'Number of primitives %d' % self.nprims)
        logger.info(self,'Maximum l in the basis %d' % self.lmax)
        logger.debug(self,'Number of primitives per center %s' % self.npc)
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
    
        title, symbols, numbers, coords, icen, ityp, oexp, \
        mo_count, mo_occ, mo_energy, mo_coeff = wfn.load_wfn(self.wfnfile)
        nprims, nmo = mo_coeff.shape

        self.natm = coords.shape[0]
        self.nelectrons = mo_occ.sum()
        self.symbols = symbols
        self.coords = coords
        self.charges = numbers
        self.nmo = nmo
        self.icen = icen
        self.ityp = ityp
        self.oexp = oexp
        self.mo_occ = mo_occ
        self.mo_coeff = mo_coeff
        self.nprims = nprims
        self.lang = numpy.empty(self.nprims, dtype=numpy.int32)
        self.npc = numpy.zeros(self.natm, dtype=numpy.int32)

        # 1) Purge repeated primitives
        purge_prims(self)

        # 2) Classifie basis set 
        get_shells(self)

        # 3) Get cutdistances for each shell
        get_shells_eps(self)

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        # 5) Test
        logger.info(self,'Testing orthogonality of orbitals')
        s = evaluator.eval_overlap(self)
        s = numpy.dot(s,self.mo_coeff)
        s = numpy.dot(self.mo_coeff.T,s)
        tools.dump_tri(self.stdout, s) 
        logger.info(self,'')

        # 5) Write info to file
        logger.info(self,'Write HDF5 file')
	mol_dic = {'natm':self.natm,
                   'symbols':self.symbols,
                   'coords':self.coords,
                   'charges':self.charges,
                   'nelectrons':self.nelectrons}
        chkfile.save(self.chkfile, 'molecule', mol_dic)
	bas_dic = {'cuttz':self.cuttz,
                   'nmo':self.nmo,
                   'nprims':self.nprims,
                   'icen':self.icen,
                   'ityp':self.ityp,
                   'oexp':self.oexp,
                   'lang':self.lang,
                   'lmax':self.lmax,
                   'npc':self.npc,
                   'mo_coeff':self.mo_coeff,
                   'mo_occ':self.mo_occ,
                   'ngroup':self.ngroup,
                   'nzexp':self.nzexp,
                   'nuexp':self.nuexp,
                   'rcutte':self.rcutte,
                   'nshells':self.nshells}
        chkfile.save(self.chkfile, 'basis', bas_dic)
        logger.timer(self,'Basis info build', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.wfn'
    mol = Mole(name)
    mol.verbose = 4
    mol.kernel()

