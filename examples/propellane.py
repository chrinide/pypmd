#!/usr/bin/env python

import numpy, sys, struct, os, avas
from pyscf import gto, scf, mcscf, symm, lib
from pyscf.tools import wfn_format

name = 'propellane'

mol = gto.Mole()
mol.atom = '''
 C1            -0.0000000000   0.0000000000   0.7780279149
 C2            -0.0000000000  -0.0000000000  -0.7780279149
 C              1.2974954650   0.0505879174   0.0000000000
 C             -0.6925581541   1.0983700753   0.0000000000
 C             -0.6049373109  -1.1489579927  -0.0000000000
 H              1.8448346224   0.9941747177  -0.0000000000
 H             -1.7833978726   1.1005862899  -0.0000000000
 H             -0.0614367498  -2.0947610077   0.0000000000
 H              1.9161504800  -0.8477829215   0.0000000000
 H             -0.2238736931   2.0833264539   0.0000000000
 H             -1.6922767869  -1.2355435324  -0.0000000000
'''
mol.basis = 'sto-3g'
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = scf.RHF(mol)
mf.max_cycle = 150
mf.kernel()

dirnow = os.path.realpath(os.path.join(__file__, '..'))
basfile = os.path.join(dirnow, 'sqzp.dat')
mol.basis = basfile
mol.build(False,False)
mf = scf.RHF(mol)
mf.max_cycle = 150
mf.chkfile = name+'.chk'
dm = mf.from_chk(name+'.chk')
#mf.__dict__.update(lib.chkfile.load(name+'.chk', 'scf'))
#dm = mf.make_rdm1()
mf.kernel(dm)
mf.analyze()

aolst1  = ['C1 2s', 'C2 2s']
aolst2  = ['C1 2p', 'C2 2p']
aolst = aolst1 + aolst2
ncas, nelecas, mo = avas.avas(mf, aolst, threshold_occ=0.1, threshold_vir=1e-5, minao='ano', ncore=5)

mc = mcscf.CASSCF(mf, ncas, nelecas)
mc.chkfile = name+'.chk'
mc.max_cycle_macro = 35
mc.max_cycle_micro = 7
#mc.fcisolver.wfnsym = 'A1g'
#mo = lib.chkfile.load(name+'.chk', 'mcscf/mo_coeff')
mc.kernel(mo)
mc.analyze()

nmo = mc.ncore + mc.ncas
rdm1, rdm2 = mc.fcisolver.make_rdm12(mc.ci, mc.ncas, mc.nelecas)
rdm1, rdm2 = mcscf.addons._make_rdm12_on_mo(rdm1, rdm2, mc.ncore, mc.ncas, nmo)

orbsym = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mc.mo_coeff[:,:nmo])
natocc, natorb = symm.eigh(-rdm1, orbsym)
for i, k in enumerate(numpy.argmax(abs(natorb), axis=0)):
    if natorb[k,i] < 0:
        natorb[:,i] *= -1
natorb = numpy.dot(mc.mo_coeff[:,:nmo], natorb)
natocc = -natocc

wfn_file = name + '.wfn'
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, natorb, mo_occ=natocc)
    wfn_format.write_coeff(f2, mol, mc.mo_coeff[:,:nmo])
    wfn_format.write_ci(f2, mc.ci, mc.ncas, mc.nelecas, ncore=mc.ncore)

