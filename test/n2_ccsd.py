#!/usr/bin/env python

import numpy, sys
from pyscf import gto, scf, cc, lib, ao2mo
from pyscf.tools import molden

name = 'n2_ccsd'

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

ncore = 2
mcc = cc.CCSD(mf)
mcc.direct = 1
mcc.diis_space = 10
mcc.frozen = ncore
mcc.conv_tol = 1e-6
mcc.conv_tol_normt = 1e-6
mcc.max_cycle = 150
mcc.kernel()

t1norm = numpy.linalg.norm(mcc.t1)
t1norm = t1norm/numpy.sqrt(mol.nelectron-ncore*2)
lib.logger.info(mcc,"* T1 norm should be les than 0.02")
lib.logger.info(mcc,"* T1 norm : %12.6f" % t1norm)

nao, nmo = mf.mo_coeff.shape
rdm1 = mcc.make_rdm1()

den_file = name + 'wfn.rdm1'
fspt = open(den_file,'w')
for i in range(nmo):
    for j in range(nmo):
        fspt.write('%i %i %.10f\n' % ((i+1), (j+1), rdm1[i,j]))
fspt.close()                    
                 
with open(name+'.mol', 'w') as f2:
    molden.header(mol, f2)
    molden.orbital_coeff(mol, f2, mf.mo_coeff[:,:nmo], occ=mf.mo_occ[:nmo])

