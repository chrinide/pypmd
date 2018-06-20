#!/usr/bin/env python

import numpy, sys
from pyscf import gto, scf
from pyscf.tools import molden

name = 'h2s_rhf'

mol = gto.Mole()
mol.basis = 'cc-pvdz'
mol.atom = '''
S      0.000000      0.000000      0.351853
H      0.000000      0.975868     -0.586476
H      0.000000     -0.975868     -0.586476
'''
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = scf.RHF(mol)
mf.kernel()

with open(name+'.mol', 'w') as f2:
    molden.header(mol, f2)
    molden.orbital_coeff(mol, f2, mf.mo_coeff[:,mf.mo_occ>0], occ=mf.mo_occ[mf.mo_occ>0])

