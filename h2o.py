#!/usr/bin/env python

import numpy, time, os, sys
from pyscf import gto, scf, lib, dft, ao2mo
from pyscf.tools import wfn_format

name = 'h2o'

mol = gto.Mole()
mol.atom = '''
O      0.000000      0.000000      0.118351
H      0.000000      0.761187     -0.469725
H      0.000000     -0.761187     -0.469725
'''
dirnow = os.path.realpath(os.path.join(__file__, '..'))
basfile = os.path.join(dirnow, 'sqzp.dat')
mol.basis = basfile
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = dft.RKS(mol)
mf.chkfile = name+'.chk'
mf.max_cycle = 150
#mf.grids.atom_grid = {'H': (20,110), 'O': (20,110)}
#mf.grids.prune = None
mf.xc = 'rpw86,pbe'
mf.kernel()

wfn_file = name + '.wfn'
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, mf.mo_coeff[:,mf.mo_occ>0], \
    mo_occ=mf.mo_occ[mf.mo_occ>0], mo_energy=mf.mo_energy[mf.mo_occ>0])

dm = mf.make_rdm1()
nao = mol.nao_nr()

coords = mf.grids.coords
weights = mf.grids.weights
ngrids = len(weights)
ao = dft.numint.eval_ao(mol, coords, deriv=1)
rho = dft.numint.eval_rho(mol, ao, dm, xctype='GGA')
lib.logger.info(mf,'Rho = %.12f' % numpy.einsum('i,i->', rho[0], weights))
ex, vx = dft.libxc.eval_xc('rPW86,', rho)[:2]
ec, vc = dft.libxc.eval_xc(',PBE', rho)[:2]
ex = numpy.einsum('i,i,i->', ex, rho[0], weights)
ec = numpy.einsum('i,i,i->', ec, rho[0], weights)

gnorm2 = numpy.zeros(ngrids)
for i in range(ngrids):
    gnorm2[i] = numpy.linalg.norm(rho[-3:,i])**2
coef_C = 0.0093
coef_B = 5.9
coef_beta = 1.0/32.0 * (3.0/(coef_B**2.0))**(3.0/4.0)
kappa_pref = coef_B * (1.5*numpy.pi)/((9.0*numpy.pi)**(1.0/6.0))
const = 4.0/3.0 * numpy.pi
vv10_e = 0.0
t = time.time()
for idx1 in range(ngrids):
    point1 = coords[idx1,:]
    rho1 = rho[0,idx1]
    weigth1 = weights[idx1]
    gamma1 = gnorm2[idx1]
    Wp1 = const*rho1
    Wg1 = coef_C * ((gamma1/(rho1*rho1))**2.0)
    W01 = numpy.sqrt(Wg1 + Wp1)
    kappa1 = rho1**(1.0/6.0)*kappa_pref
    #
    R =  (point1[0]-coords[:,0])**2
    R += (point1[1]-coords[:,1])**2
    R += (point1[2]-coords[:,2])**2
    Wp2 = const*rho[0]
    Wg2 = coef_C * ((gnorm2/(rho[0]*rho[0]))**2.0)
    W02 = numpy.sqrt(Wg2 + Wp2)
    kappa2 = rho[0]**(1.0/6.0)*kappa_pref
    g = W01*R + kappa1
    gp = W02*R + kappa2
    kernel12 = -1.5*weights*rho[0]/(g*gp*(g+gp))
    # Energy 
    kernel = kernel12.sum()
    vv10_e += weigth1*rho1*(coef_beta + 0.5*kernel)
lib.logger.info(mf,'VV10 = %.12f' % vv10_e)
lib.logger.info(mf,'Total time taken VV10: %.3f seconds' % (time.time()-t))

s = mol.intor('int1e_ovlp')
t = mol.intor('int1e_kin')
v = mol.intor('int1e_nuc')
eri_ao = ao2mo.restore(1,mf._eri,nao)
eri_ao = eri_ao.reshape(nao,nao,nao,nao)

enuc = mol.energy_nuc() 
ekin = numpy.einsum('ij,ji->',t,dm)
pop = numpy.einsum('ij,ji->',s,dm)
elnuce = numpy.einsum('ij,ji->',v,dm)
lib.logger.info(mf,'Population : %12.6f' % pop)
lib.logger.info(mf,'Kinetic energy : %12.6f' % ekin)
lib.logger.info(mf,'Nuclear Atraction energy : %12.6f' % elnuce)
lib.logger.info(mf,'Nuclear Repulsion energy : %12.6f' % enuc)
bie1 = numpy.einsum('ijkl,ij,kl->',eri_ao,dm,dm)*0.5 # J
bie2 = numpy.einsum('ijkl,il,jk->',eri_ao,dm,dm)*0.25 # XC
pairs1 = numpy.einsum('ij,kl,ij,kl->',dm,dm,s,s) # J
pairs2 = numpy.einsum('ij,kl,li,kj->',dm,dm,s,s)*0.5 # XC
pairs = (pairs1 - pairs2)
lib.logger.info(mf,'Coulomb Pairs : %12.6f' % (pairs1))
lib.logger.info(mf,'XC Pairs : %12.6f' % (pairs2))
lib.logger.info(mf,'Pairs : %12.6f' % pairs)
lib.logger.info(mf,'J energy : %12.6f' % bie1)
lib.logger.info(mf,'XC energy : %12.6f' % -bie2)
lib.logger.info(mf,'EE energy : %12.6f' % (bie1-bie2))
etot = enuc + ekin + elnuce + bie1 - bie2
lib.logger.info(mf,'HF Total energy : %12.6f' % etot)
lib.logger.info(mf,'Ex : %12.6f' % ex)
lib.logger.info(mf,'Ec : %12.6f' % ec)
lib.logger.info(mf,'Exc : %12.6f' % (ex+ec))
etot = enuc + ekin + elnuce + bie1 + ex + ec
lib.logger.info(mf,'DFT Total energy : %12.6f' % etot)

