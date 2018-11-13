#!/usr/bin/env python

import os, numpy
from pyscf import gto, scf
from pyscf.tools import wfn_format
from pyscf.dft import numint

name = 'h2o'

mol = gto.Mole()
mol.atom = '''
  O  0.0000  0.0000  0.1173
  H  0.0000  0.7572 -0.4692
  H  0.0000 -0.7572 -0.4692
'''
dirnow = os.path.realpath(os.path.join(__file__, '..'))
basfile = os.path.join(dirnow, 'sqzp.dat')
mol.basis = basfile
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = scf.RHF(mol)
mf.max_cycle = 150
mf.chkfile = name+'.chk'
#mf.__dict__.update(lib.chkfile.load(name+'.chk', 'scf'))
#dm = mf.make_rdm1()
mf.kernel()
dm = mf.make_rdm1()

coeff = mf.mo_coeff[:,mf.mo_occ>0]
occ = mf.mo_occ[mf.mo_occ>0]
energy = mf.mo_energy[mf.mo_occ>0]
wfn_file = name + '.wfn'
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, coeff, mo_occ=occ, mo_energy=energy)

coords = numpy.random.random((5,3))*3.0
print "Coordinates and density values"
ao_value = numint.eval_ao(mol, coords, deriv=1)
# The first row of rho is electron density, the rest three rows are electron
# density gradients which are needed for GGA functional
rho = numint.eval_rho(mol, ao_value, dm, xctype='GGA')
print coords 
print rho
