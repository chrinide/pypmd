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
#

import sys
import h5py
import time
import numpy
import ctypes

import misc
import param
import logger
import chkfile

# For code compatiblity in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

libfapi = misc.load_library('libfapi')
libcapi = misc.load_library('libcapi')

# ~= (L+1)**2/3
LEBEDEV_ORDER = {
      0:    1,
      3:    6,
      5:   14,
      7:   26,
      9:   38,
     11:   50,
     13:   74,
     15:   86,
     17:  110,
     19:  146,
     21:  170,
     23:  194,
     25:  230,
     27:  266,
     29:  302,
     31:  350,
     35:  434,
     41:  590,
     47:  770,
     53:  974,
     59: 1202,
     65: 1454,
     71: 1730,
     77: 2030,
     83: 2354,
     89: 2702,
     95: 3074,
    101: 3470,
    107: 3890,
    113: 4334,
    119: 4802,
    125: 5294,
    131: 5810
}

LEBEDEV_NGRID = numpy.asarray((
    1   , 6   , 14  , 26  , 38  , 50  , 74  , 86  , 110 , 146 ,
    170 , 194 , 230 , 266 , 302 , 350 , 434 , 590 , 770 , 974 ,
    1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334,
    4802, 5294, 5810))

#########################
# JCP 41 3199 (1964).
BRAGG = 1.0/param.BOHR * numpy.array((0,  # Ghost atom
        0.35,                                     1.40,             # 1s
        1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 1.50,             # 2s2p
        1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.80,             # 3s3p
        2.20, 1.80,                                                 # 4s
        1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, # 3d
                    1.30, 1.25, 1.15, 1.15, 1.15, 1.90,             # 4p
        2.35, 2.00,                                                 # 5s
        1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, # 4d
                    1.55, 1.45, 1.45, 1.40, 1.40, 2.10,             # 5p
        2.60, 2.15,                                                 # 6s
        1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85,                   # La, Ce-Eu
        1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,             # Gd, Tb-Lu
              1.55, 1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35, 1.50, # 5d
                    1.90, 1.80, 1.60, 1.90, 1.45, 2.10,             # 6p
        1.80, 2.15,                                                 # 7s
        1.95, 1.80, 1.80, 1.75, 1.75, 1.75, 1.75,
        1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
        1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
                    1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
        1.75, 1.75,
        1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75))

# from Gerald Knizia's CtDftGrid, which is based on
#       http://en.wikipedia.org/wiki/Covalent_radius
# and
#       Beatriz Cordero, Veronica Gomez, Ana E. Platero-Prats, Marc Reves,
#       Jorge Echeverria, Eduard Cremades, Flavia Barragan and Santiago
#       Alvarez.  Covalent radii revisited. Dalton Trans., 2008, 2832-2838,
#       doi:10.1039/b801115j
COVALENT = 1.0/param.BOHR * numpy.array((0,  # Ghost atom
        0.31,                                     0.28,             # 1s
        1.28, 0.96, 0.84, 0.73, 0.71, 0.66, 0.57, 0.58,             # 2s2p
        1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,             # 3s3p
        2.03, 1.76,                                                 # 4s
        1.70, 1.60, 1.53, 1.39, 1.50, 1.42, 1.38, 1.24, 1.32, 1.22, # 3d
                    1.22, 1.20, 1.19, 1.20, 1.20, 1.16,             # 4p
        2.20, 1.95,                                                 # 5s
        1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, # 4d
                    1.42, 1.39, 1.39, 1.38, 1.39, 1.40,             # 5p
        2.44, 2.15,                                                 # 6s
        2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98,                   # La, Ce-Eu
        1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,             # Gd, Tb-Lu
              1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, # 5d
                    1.45, 1.46, 1.48, 1.40, 1.50, 1.50,             # 6p
        2.60, 2.21,                                                 # 7s
        2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69))

def stratmann(g):
    '''Stratmann, Scuseria, Frisch. CPL, 257, 213 (1996)'''
    a = 0.64  # for eq. 14
    g = numpy.asarray(g)
    ma = g/a
    ma2 = ma * ma
    g1 = numpy.asarray((1.0/16.0)*(ma*(35.0 + ma2*(-35.0 + ma2*(21.0 - 5.0*ma2)))))
    g1[g<=-a] = -1.0
    g1[g>= a] =  1.0
    return g1

def original_becke(g):
    '''Becke, JCP, 88, 2547 (1988)'''
    g = (3.0 - g**2) * g*0.5
    g = (3.0 - g**2) * g*0.5
    g = (3.0 - g**2) * g*0.5
    return g

# Gauss-Chebyshev of the first kind,  and the transformed interval [0,\infty)
def becke(n, charge, *args, **kwargs):
    '''Becke, JCP, 88, 2547 (1988)'''
    if charge == 1:
        rm = BRAGG[charge]
    else:
        rm = BRAGG[charge] * 0.5
    t, w = numpy.polynomial.chebyshev.chebgauss(n)
    r = (1.0+t)/(1.0-t) * rm
    w *= 2.0/(1.0-t)**2 * rm
    return r, w

# scale rad and rad_weight if necessary
# gauss-legendre
def delley(n, *args, **kwargs):
    '''B. Delley radial grids. Ref. JCP, 104, 9848 log2 algorithm'''
    r = numpy.empty(n)
    dr = numpy.empty(n)
    r_outer = 12.0
    step = 1.0 / (n+1)
    rfac = r_outer / numpy.log(1.0-(n*step)**2)
    for i in range(1, n+1):
        xi = rfac * numpy.log(1.0-(i*step)**2)
        r[i-1] = xi
        dri = rfac * (-2.0*i*(step)**2) / ((1.0-(i*step)**2)) # d xi / dr
        dr[i-1] = dri
    return r, dr
gauss_legendre = delley

def mura_knowles(n, charge=None, *args, **kwargs):
    '''Mura-Knowles (JCP, 104, 9848) log3 quadrature radial grids'''
    r = numpy.empty(n)
    dr = numpy.empty(n)
# 7 for Li, Be, Na, Mg, K, Ca, otherwise 5
    if charge in (3, 4, 11, 12, 19, 20):
        far = 7
    else:
        far = 5.2
    for i in range(n):
        x = (i+0.5) / n
        r[i] = -far * numpy.log(1.0-x**3)
        dr[i] = far * 3.0*x*x/((1.0-x**3)*n)
    return r, dr

# Gauss-Chebyshev of the second kind,  and the transformed interval [0,\infty)
# Ref  Matthias Krack and Andreas M. Koster,  J. Chem. Phys. 108 (1998), 3226
def gauss_chebyshev(n, *args, **kwargs):
    '''Gauss-Chebyshev (JCP, 108, 3226) radial grids'''
    ln2 = 1.0 / numpy.log(2.0)
    fac = 16.0/3.0 / (n+1)
    x1 = numpy.arange(1,n+1) * numpy.pi / (n+1)
    xi = ((n-1-numpy.arange(n)*2) / (n+1.0) +
          (1.0+2.0/3.0*numpy.sin(x1)**2) * numpy.sin(2.0*x1) / numpy.pi)
    xi = (xi - xi[::-1])/2.0
    r = 1.0 - numpy.log(1.0+xi) * ln2
    dr = fac * numpy.sin(x1)**4 * ln2/(1.0+xi)
    return r, dr

def treutler_ahlrichs(n, *args, **kwargs):
    '''
    Treutler-Ahlrichs (JCP 102, 346 (M4)) radial grids
    '''
    r = numpy.empty(n)
    dr = numpy.empty(n)
    step = numpy.pi / (n+1)
    ln2 = 1.0 / numpy.log(2.0)
    for i in range(n):
        x = numpy.cos((i+1)*step)
        r[i] = -ln2*(1.0+x)**0.6 * numpy.log((1.0-x)/2.0)
        dr[i] = step * numpy.sin((i+1.0)*step) \
                * ln2*(1.0+x)**0.6 *(-0.6/(1.0+x)*numpy.log((1.0-x)/2.0)+1.0/(1.0-x))
    return r[::-1], dr[::-1]
treutler = treutler_ahlrichs

def becke_atomic_radii_adjust(self, atomic_radii):
    '''Becke atomic radii adjust function'''
# Becke atomic size adjustment.  J. Chem. Phys. 88, 2547
# i > j
# fac(i,j) = \frac{1}{4} ( \frac{ra(j)}{ra(i)} - \frac{ra(i)}{ra(j)}
# fac(j,i) = -fac(i,j)

    charges = self.charges
    rad = atomic_radii[charges] + 1e-200
    rr = rad.reshape(-1,1) * (1.0/rad)
    a = 0.25 * (rr.T - rr)
    a[a<-0.5] = -0.5
    a[a>0.5] = 0.5
    #:return lambda i,j,g: g + a[i,j]*(1-g**2)
    def fadjust(i, j, g):
        g1 = g**2
        g1 -= 1.0
        g1 *= -a[i,j]
        g1 += g
        return g1
    return fadjust

def treutler_atomic_radii_adjust(self, atomic_radii):
    '''Treutler atomic radii adjust function: JCP, 102, 346'''
# JCP, 102, 346
# i > j
# fac(i,j) = \frac{1}{4} ( \frac{ra(j)}{ra(i)} - \frac{ra(i)}{ra(j)}
# fac(j,i) = -fac(i,j)
    charges = self.charges
    rad = numpy.sqrt(atomic_radii[charges]) + 1e-200
    rr = rad.reshape(-1,1) * (1.0/rad)
    a = 0.25 * (rr.T - rr)
    a[a<-0.5] = -0.5
    a[a>0.5] = 0.5
    #:return lambda i,j,g: g + a[i,j]*(1-g**2)
    def fadjust(i, j, g):
        g1 = g**2
        g1 -= 1.0
        g1 *= -a[i,j]
        g1 += g
        return g1
    return fadjust

def _default_rad(nuc, level=3):
    '''Number of radial grids '''
    tab   = numpy.array( (2 , 10, 18, 36, 54, 86, 118))
    #           Period    1   2   3   4   5   6   7         # level
    grids = numpy.array((( 10, 15, 20, 30, 35, 40, 50),     # 0
                         ( 30, 40, 50, 60, 65, 70, 75),     # 1
                         ( 40, 60, 65, 75, 80, 85, 90),     # 2
                         ( 50, 75, 80, 90, 95,100,105),     # 3
                         ( 60, 90, 95,105,110,115,120),     # 4
                         ( 70,105,110,120,125,130,135),     # 5
                         ( 80,120,125,135,140,145,150),     # 6
                         ( 90,135,140,150,155,160,165),     # 7
                         (100,150,155,165,170,175,180),     # 8
                         (200,200,200,200,200,200,200),))   # 9
    period = (nuc > tab).sum()
    return grids[level,period]

def _default_ang(nuc, level=3):
    '''Order of angular grids. See LEBEDEV_ORDER for the mapping of
    the order and the number of angular grids'''
    tab   = numpy.array( (2 , 10, 18, 36, 54, 86, 118))
    #           Period    1   2   3   4   5   6   7         # level
    order = numpy.array(((11, 15, 17, 17, 17, 17, 17 ),     # 0
                         (17, 23, 23, 23, 23, 23, 23 ),     # 1
                         (23, 29, 29, 29, 29, 29, 29 ),     # 2
                         (29, 29, 35, 35, 35, 35, 35 ),     # 3
                         (35, 41, 41, 41, 41, 41, 41 ),     # 4
                         (41, 47, 47, 47, 47, 47, 47 ),     # 5
                         (47, 53, 53, 53, 53, 53, 53 ),     # 6
                         (53, 59, 59, 59, 59, 59, 59 ),     # 7
                         (59, 59, 59, 59, 59, 59, 59 ),     # 8
                         (65, 65, 65, 65, 65, 65, 65 ),))   # 9
    period = (nuc > tab).sum()
    return LEBEDEV_ORDER[order[level,period]]

def lebgrid(npang):

    if npang not in LEBEDEV_NGRID:
        raise ValueError('Lebgrid unsupported angular grid %d' % npang)
    else:
        ct = numpy.zeros((npang), order='F')
        st = numpy.zeros((npang), order='F') 
        cp = numpy.zeros((npang), order='F') 
        sp = numpy.zeros((npang), order='F') 
        aw = numpy.zeros((npang), order='F') 
        agrids = numpy.zeros((npang,5))
        libfapi.lebgrid(ct.ctypes.data_as(ctypes.c_void_p), 
                        st.ctypes.data_as(ctypes.c_void_p), 
                        cp.ctypes.data_as(ctypes.c_void_p), 
                        sp.ctypes.data_as(ctypes.c_void_p), 
                        aw.ctypes.data_as(ctypes.c_void_p),  
                        ctypes.c_int(npang)) 
        agrids[:,0] = ct
        agrids[:,1] = st
        agrids[:,2] = cp
        agrids[:,3] = sp
        agrids[:,4] = aw

    return agrids

def prange(start, end, step):
    for i in range(start, end, step):
        yield i, min(i+step, end)

def gen_atomic_grids(self, atom_grid={}, radi_method=gauss_chebyshev, level=3, **kwargs):
    '''Generate number of radial grids and angular grids for the given molecule.

    Returns:
        A dict, with the atom symbol for the dict key.  For each atom type,
        the dict value has two items: one is the meshgrid coordinates wrt the
        atom center; the second is the volume of that grid.
    '''
    if isinstance(atom_grid, (list, tuple)):
        atom_grid = dict([(self.symbols[ia], atom_grid)
                          for ia in range(self.natm)])
    atom_grids_tab = {}
    for ia in range(self.natm):
        symb = self.symbols[ia]
        if symb not in atom_grids_tab:
            #chg = self.charges[symb] #Better way in case of ECP of rare symbol
            chg = self.charges[ia]
            if symb in atom_grid:
                n_rad, n_ang = atom_grid[symb]
                if n_ang not in LEBEDEV_NGRID:
                    if n_ang in LEBEDEV_ORDER:
                        logger.warn(self, 'n_ang %d for atom %d %s is not '
                                    'the supported Lebedev angular grids. '
                                    'Set n_ang to %d', n_ang, ia, symb,
                                    LEBEDEV_ORDER[n_ang])
                        n_ang = LEBEDEV_ORDER[n_ang]
                    else:
                        raise ValueError('Unsupported angular grids %d' % n_ang)
            else:
                n_rad = _default_rad(chg, level)
                n_ang = _default_ang(chg, level)
            rad, dr = radi_method(n_rad, chg, ia, **kwargs)
            rad_weight = 4.0*numpy.pi * rad**2 * dr
            angs = [n_ang] * n_rad
            logger.debug(self,'atom %s rad-grids = %d, ang-grids = %s', symb, n_rad, angs)

            angs = numpy.array(angs)
            coords = []
            vol = []
            for n in sorted(set(angs)):
                grid = numpy.empty((n,4))
                libcapi.MakeAngularGrid(grid.ctypes.data_as(ctypes.c_void_p),
                                        ctypes.c_int(n))
                idx = numpy.where(angs==n)[0]
                for i0, i1 in prange(0, len(idx), 12):  # 12 radi-grids as a group
                    coords.append(numpy.einsum('i,jk->jik',rad[idx[i0:i1]],
                                               grid[:,:3]).reshape(-1,3))
                    vol.append(numpy.einsum('i,j->ji', rad_weight[idx[i0:i1]],
                                            grid[:,3]).ravel())
            atom_grids_tab[symb] = (numpy.vstack(coords), numpy.hstack(vol))
    return atom_grids_tab

def inter_distance(self):
    '''
    Inter-particle distance array
    '''
    coords = self.xyz
    rr = numpy.linalg.norm(coords.reshape(-1,1,3) - coords, axis=2)
    rr[numpy.diag_indices_from(rr)] = 0
    return rr

def gen_partition(self, atom_grids_tab,
                  radii_adjust=None, atomic_radii=BRAGG,
                  becke_scheme=original_becke):
    '''Generate the mesh grid coordinates and weights for DFT numerical integration.
    We can change radii_adjust, becke_scheme functions to generate different meshgrid.

    Returns:
        grid_coord and grid_weight arrays.  grid_coord array has shape (N,3);
        weight 1D array has N elements.
    '''
    if callable(radii_adjust) and atomic_radii is not None:
        f_radii_adjust = radii_adjust(self, atomic_radii)
    else:
        f_radii_adjust = None
    atm_coords = numpy.asarray(self.xyz , order='C')
    atm_dist = inter_distance(self)
    def gen_grid_partition(coords):
        ngrids = coords.shape[0]
        grid_dist = numpy.empty((self.natm,ngrids))
        for ia in range(self.natm):
            dc = coords - atm_coords[ia]
            grid_dist[ia] = numpy.sqrt(numpy.einsum('ij,ij->i',dc,dc))
        pbecke = numpy.ones((self.natm,ngrids))
        for i in range(self.natm):
            for j in range(i):
                g = 1.0/atm_dist[i,j] * (grid_dist[i]-grid_dist[j])
                if f_radii_adjust is not None:
                    g = f_radii_adjust(i, j, g)
                g = becke_scheme(g)
                pbecke[i] *= 0.5 * (1.0-g)
                pbecke[j] *= 0.5 * (1.0+g)
        return pbecke

    coords_all = []
    weights_all = []
    for ia in range(self.natm):
        coords, vol = atom_grids_tab[self.symbols[ia]]
        coords = coords + atm_coords[ia]
        pbecke = gen_grid_partition(coords)
        weights = vol * pbecke[ia] * (1.0/pbecke.sum(axis=0))
        coords_all.append(coords)
        weights_all.append(weights)
    return numpy.vstack(coords_all), numpy.hstack(weights_all)


class Grids(object):
    '''DFT mesh grids

    Attributes for Grids:
        level : int (0 - 9)
            big number for large mesh grids, default is 3

        atomic_radii : 1D array
            | BRAGG  (default)
            | COVALENT
            | None : to switch off atomic radii adjustment

        radii_adjust : function(self, atomic_radii) => (function(atom_id, atom_id, g) => array_like_g)
            Function to adjust atomic radii, can be one of
            | treutler_atomic_radii_adjust
            | becke_atomic_radii_adjust
            | None : to switch off atomic radii adjustment

        radi_method : function(n) => (rad_grids, rad_weights)
            scheme for radial grids, can be one of
            | treutler  (default)
            | delley
            | mura_knowles
            | gauss_chebyshev

        becke_scheme : function(v) => array_like_v
            weight partition function, can be one of
            | original_becke  (default)
            | stratmann

        atom_grid : dict
            Set (radial, angular) grids for particular atoms.
            Eg, grids.atom_grid = {'H': (20,110)} will generate 20 radial
            grids and 110 angular grids for H atom.

        level : int
            To control the number of radial and angular grids.  The default
            level 3 corresponds to
            (50,302) for H, He;
            (75,302) for second row;
            (80~105,434) for rest.

        '''
    def __init__(self, datafile):
        self.verbose = logger.NOTE
        self.chkfile = datafile
        self.stdout = sys.stdout
        self.max_memory = param.MAX_MEMORY
        self.chkfile = datafile
        self.scratch = param.TMPDIR 
        self.level = 3
        self.atom_grid = {}
        self.atomic_radii = BRAGG
        self.radii_adjust = treutler_atomic_radii_adjust
        self.radi_method =  treutler
        self.becke_scheme = original_becke
##################################################
# don't modify the following attributes, they are not input options
        self.xyz = None
        self.charges = None
        self.symbols = None
        self.natm = None
        self.coords  = None
        self.weights = None
        self._keys = set(self.__dict__.keys())

    def __setattr__(self, key, val):
        if key in ('atom_grid', 'atomic_radii', 'radii_adjust', 'radi_method',
                   'becke_scheme', 'level'):
            self.coords = None
            self.weights = None
            self.non0tab = None
        super(Grids, self).__setattr__(key, val)

    def dump_flags(self):

        if self.verbose < logger.INFO:
            return self

        logger.info(self,'')
        logger.info(self,'******** %s flags ********', self.__class__)
        logger.info(self,'* General Info')
        logger.info(self,'Verbose level %d' % self.verbose)
        logger.info(self,'Scratch dir %s' % self.scratch)
        logger.info(self,'Input h5 data file %s' % self.chkfile)

        logger.info(self,'* Molecular Info')
        logger.info(self,'Num atoms %d' % self.natm)
        logger.info(self,'Atom Coordinates (Bohr)')
        for i in range(self.natm):
            logger.info(self,'Nuclei %d %s with charge %d position : %.6f  %.6f  %.6f', i, 
                        self.symbols[i], self.charges[i], *self.xyz[i])

        logger.info(self,'* Grid Info')
        logger.info(self,'radial grids: %s', self.radi_method)
        logger.info(self,'becke partition: %s', self.becke_scheme)
        logger.info(self,'grids dens level: %d', self.level)
        if self.radii_adjust is not None:
            logger.info(self,'atomic radii adjust function: %s',
                        self.radii_adjust)
            logger.debug(self,'atomic_radii : %s', self.atomic_radii)
        if self.atom_grid:
            logger.info(self,'User specified grid scheme %s', str(self.atom_grid))
        logger.info(self,'')

        return self

    def build(self, **kwargs):

        t0 = time.clock()
        logger.TIMER_LEVEL = 3

        # 1) Read info
        logger.info(self,'Reading HDF5 file')
        with h5py.File(self.chkfile) as f:
            self.natm = f['molecule/natm'].value
            self.xyz = f['molecule/coords'].value
            self.charges = f['molecule/charges'].value
            self.symbols = f['molecule/symbols'].value

        # 2) Dump info
        if self.verbose > logger.NOTE:
            self.dump_flags()

        # 3) Gen grids
        atom_grids_tab = gen_atomic_grids(self,self.atom_grid,
                                               self.radi_method,
                                               self.level, **kwargs)
        self.coords, self.weights = \
                gen_partition(self,atom_grids_tab,
                                   self.radii_adjust, self.atomic_radii,
                                   self.becke_scheme)
        # 4) Save
        logger.info(self,'Finish with grids')
        logger.info(self,'Tot grids %d', len(self.weights))
        logger.info(self,'Write HDF5 grid file')
        grid_dic = {'coords':self.coords,
                    'weights':self.weights}
        chkfile.save(self.chkfile, 'becke', grid_dic)
        logger.info(self,'Grid saved')
        logger.timer(self,'Becke grid build', t0)
        logger.info(self,'')

        return self

    kernel = build

if __name__ == '__main__':

    npang = 5810
    agrid = lebgrid(npang)
    with open('agrid.txt', 'w') as f2:
        f2.write('# Point 1 2 3 4 weight\n')
        for i in range(npang):
            f2.write('%d   %.6f  %.6f  %.6f  %.6f  %.6f\n' % \
            ((i+1), agrid[i,0], agrid[i,1], agrid[i,2], agrid[i,3], agrid[i,4]))

    name = 'h2o.wfn.h5'
    g = Grids(name)
    g.verbose = 4
    g.build()
    print(g.coords.shape)
    with open('bgrid.txt', 'w') as f2:
        f2.write('# Point weight\n')
        for i in range(g.weights.shape[0]):
            f2.write('%d   %.6f  %.6f  %.6f  %.6f\n' % \
            ((i+1), g.coords[i,0], g.coords[i,1], g.coords[i,2], g.weights[i]))
