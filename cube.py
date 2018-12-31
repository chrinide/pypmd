#!/usr/bin/env python
# Copyright 2014-2018 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#         Peter Koval <koval.peter@gmail.com>
#         Paul J. Robinson <pjrobinson@ucla.edu>
#

'''
Gaussian cube file format.  Reference:
http://paulbourke.net/dataformats/cube/
http://gaussian.com/cubegen/

The output cube file has the following format

Comment line
Comment line
N_atom Ox Oy Oz         # number of atoms, followed by the coordinates of the origin
N1 vx1 vy1 vz1          # number of grids along each axis, followed by the step size in x/y/z direction.
N2 vx2 vy2 vz2          # ...
N3 vx3 vy3 vz3          # ...
Atom1 Z1 x y z          # Atomic number, charge, and coordinates of the atom
...                     # ...
AtomN ZN x y z          # ...
Data on grids           # (N1*N2) lines of records, each line has N3 elements
'''

import sys
import time
import h5py
import numpy

import param
import data
import numint
import logger

RESOLUTION = None
BOX_MARGIN = 3.0

def cartesian_prod(arrays, out=None):
    '''
    Generate a cartesian product of input arrays.
    http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

    Args:
        arrays : list of array-like
            1-D arrays to form the cartesian product of.
        out : ndarray
            Array to place the cartesian product in.

    Returns:
        out : ndarray
            2-D array of shape (M, len(arrays)) containing cartesian products
            formed of input arrays.

    Examples:

    >>> cartesian_prod(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    '''
    arrays = [numpy.asarray(x) for x in arrays]
    dtype = numpy.result_type(*arrays)
    nd = len(arrays)
    dims = [nd] + [len(x) for x in arrays]

    if out is None:
        out = numpy.empty(dims, dtype)
    else:
        out = numpy.ndarray(dims, dtype, buffer=out)
    tout = out.reshape(dims)

    shape = [-1] + [1] * nd
    for i, arr in enumerate(arrays):
        tout[i] = arr.reshape(shape[:nd-i])

    return tout.reshape(nd,-1).T

def prange(start, end, step):
    for i in range(start, end, step):
        yield i, min(i+step, end)

def density(self):
    """Calculates electron density and write out in cube format.
    """
    # Compute density on the cube grid
    coords = self.get_coords()
    ngrids = self.get_ngrids()
    blksize = min(8000, ngrids)
    rho = numpy.empty(ngrids)
    for ip0, ip1 in prange(0, ngrids, blksize):
        rho[ip0:ip1] = numint.eval_rho(self,coords[ip0:ip1])
    rho = rho.reshape(self.nx,self.ny,self.nz)
    # Write out density to the .cube file
    outfile = self.chkfile+'.cube'
    self.write(rho, outfile, comment='Electron density in real space (e/Bohr^3)')

class Cube(object):
    '''  Read-write of the Gaussian CUBE files  '''
    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
        self.nx = 80
        self.ny = 80
        self.nz = 80
        self.resolution = None
        self.margin = BOX_MARGIN
        self.gpu = False
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
        self.xs = None
        self.ys = None
        self.zs = None
        self.box = None
        self.boxorig = None
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

        logger.info(self,'* Cube Info')
        logger.info(self,'Grid resolution %s', self.resolution)
        logger.info(self,'Number of nx,ny,nz %d %d %d' % (self.nx,self.ny,self.nz))
        logger.info(self,'Margin %d' % self.margin)
        logger.info(self,'Box origin %.6f %.6f %.6f', *self.boxorig)
        #self.box = None
        logger.info(self,'')

        return self

    def build(self):

        t0 = time.clock()
        logger.TIMER_LEVEL = 3
    
        # 1) Read info
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

        # 2) Build cube data
        coords = self.coords
        margin = self.margin
        resolution = self.resolution
        self.box = numpy.max(coords,axis=0) - numpy.min(coords,axis=0) + margin*2
        self.boxorig = numpy.min(coords,axis=0) - margin
        nx = self.nx
        ny = self.ny
        nz = self.nz
        if resolution is not None:
            nx, ny, nz = numpy.ceil(self.box / resolution).astype(int)
            self.nx = nx
            self.ny = ny
            self.nz = nz
        # .../(nx-1) to get symmetric mesh
        # see also the discussion on https://github.com/sunqm/pyscf/issues/154
        self.xs = numpy.arange(nx) * (self.box[0] / (nx - 1))
        self.ys = numpy.arange(ny) * (self.box[1] / (ny - 1))
        self.zs = numpy.arange(nz) * (self.box[2] / (nz - 1))

        # Dump info
        if self.verbose > logger.NOTE:
            self.dump_input()

        density(self)

        logger.timer(self,'Cube done', t0)
        logger.info(self,'')

        return self

    kernel = build

    def get_coords(self) :
        """  Result: set of coordinates to compute a field which is to be stored
        in the file.
        """
        coords = cartesian_prod([self.xs,self.ys,self.zs])
        coords = numpy.asarray(coords, order='C') - (-self.boxorig)
        return coords

    def get_ngrids(self):
        return self.nx * self.ny * self.nz

    def get_volume_element(self):
        return (self.xs[1]-self.xs[0])*(self.ys[1]-self.ys[0])*(self.zs[1]-self.zs[0])

    def write(self, field, fname, comment=None):
        """  Result: .cube file with the field in the file fname.  """
        assert(field.ndim == 3)
        assert(field.shape == (self.nx, self.ny, self.nz))
        if comment is None:
            comment = 'Generic field? Supply the optional argument "comment" to define this line'

        coords = self.coords
        with open(fname, 'w') as f:
            f.write(comment+'\n')
            f.write('PyPMD Date: %s\n' % (time.ctime()))
            f.write('%5d' % self.natm)
            f.write('%12.6f%12.6f%12.6f\n' % tuple(self.boxorig.tolist()))
            f.write('%5d%12.6f%12.6f%12.6f\n' % (self.nx, self.xs[1], 0, 0))
            f.write('%5d%12.6f%12.6f%12.6f\n' % (self.ny, 0, self.ys[1], 0))
            f.write('%5d%12.6f%12.6f%12.6f\n' % (self.nz, 0, 0, self.zs[1]))
            for ia in range(self.natm):
                chg = self.charges[ia]
                f.write('%5d%12.6f'% (chg, chg))
                f.write('%12.6f%12.6f%12.6f\n' % tuple(coords[ia]))

            for ix in range(self.nx):
                for iy in range(self.ny):
                    for iz0, iz1 in prange(0, self.nz, 6):
                        fmt = '%13.5E' * (iz1-iz0) + '\n'
                        f.write(fmt % tuple(field[ix,iy,iz0:iz1].tolist()))

#del(RESOLUTION, BOX_MARGIN)

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    cube = Cube(name)
    cube.verbose = 4
    cube.kernel()

