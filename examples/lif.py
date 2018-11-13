#!/usr/bin/env python

from pyscf import gto, scf, lib
from pyscf.tools import wfn_format

name = 'lif'

mol = gto.Mole()
mol.atom = '''
  F     0.0000  0.0000  0.0000
  Li    0.0000  0.0000  1.5639
'''
mol.basis = 'aug-cc-pvdz'
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = scf.RHF(mol)
mf.max_cycle = 150
mf.chkfile = name+'.chk'
mf.kernel()

wfn_file = name + '.wfn'
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, mf.mo_coeff[:,mf.mo_occ>0], \
    mo_occ=mf.mo_occ[mf.mo_occ>0])

