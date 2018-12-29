#!/usr/bin/env python

import numpy, time
from pyscf import gto, scf, lib, dft
from pyscf.tools import wfn_format

name = 'h2o'

mol = gto.Mole()
mol.atom = '''
O      0.000000      0.000000      0.118351
H      0.000000      0.761187     -0.469725
H      0.000000     -0.761187     -0.469725
'''
mol.basis = 'sto-6g'
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = dft.RKS(mol)
mf.chkfile = name+'.chk'
mf.max_cycle = 150
mf.grids.atom_grid = {'H': (20,110), 'O': (20,110)}
mf.grids.prune = None
mf.xc = 'rpw86,pbe'
mf.kernel()

wfn_file = name + '.wfn'
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, mf.mo_coeff[:,mf.mo_occ>0], \
    mo_occ=mf.mo_occ[mf.mo_occ>0])

