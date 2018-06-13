#!/usr/bin/env python

import numpy, sys
from pyscf import gto, scf
from pyscf.tools import molden

name = 'n2_rhf'

mol = gto.Mole()
mol.basis = 'cc-pvdz'
mol.atom = '''
N  0.0000  0.0000  0.5488
N  0.0000  0.0000 -0.5488
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

