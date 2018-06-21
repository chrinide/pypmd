#!/usr/bin/env python

import numpy, sys
from pyscf import gto, scf
from pyscf.tools import molden

name = 'h2o_h2o_rhf'

mol = gto.Mole()
mol.basis = 'cc-pvdz'
mol.atom = '''
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561
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

