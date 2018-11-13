#!/usr/bin/python

import numpy
import sys
import time

from pypmd.lib import data

def rhograd(point,self):

    fun = numpy.zeros((3))
    fun1 = numpy.zeros((3))
    gun = numpy.zeros((self.nmo))
    gun1 = numpy.zeros((self.nmo,3))
    xcoor = numpy.zeros((3))

    for i in range(self.nprim):
        ic = self.icen[i]-1
        itip = self.ityp[i]-1
        it = data.nlm[itip,:]
        ori = -self.oexp[i]
        dp2 = ori+ori
        xcoor = point - self.coords[ic,:]
        dis2 = numpy.einsum('i,i->', xcoor, xcoor)
        aexp = numpy.exp(ori*dis2)
        for j in range(3):
            n = it[j]
            x = xcoor[j]
            if (n == 0):
                dp2x = dp2*x
                fun1[j] = dp2x
                fun[j] = 1.0
            elif (n == 1):
                x2 = x*x
                dp2x2 = dp2*x2
                fun1[j] = 1.0+dp2x2
                fun[j] = x
            elif (n == 2):
                x2 = x*x
                dp2x2 = dp2*x2
                fun1[j] = x*(2.0+dp2x2)
                fun[j] = x2
            elif (n == 3):
                x2 = x*x
                dp2x2 = dp2*x2
                fun1[j] = x2*(3.0+dp2x2)
                fun[j] = x*x2
            elif (n == 4):
                x2 = x*x
                dp2x2 = dp2*x2
                fun1[j] = x2*x*(4.0+dp2x2)
                fun[j] = x2*x2
            elif (n == 5):
                x2 = x*x
                dp2x2 = dp2*x2
                fun1[j] = x2*x2*(5.0+dp2x2)
                fun[j] = x2*x2*x

        f23 = fun[1]*fun[2]*aexp
        f13 = fun[0]*fun[2]*aexp
        f12 = fun[0]*fun[1]*aexp
        gun += self.mo_coeff[i,:]*fun[0]*f23
        gun1[:,0] += self.mo_coeff[i,:]*fun1[0]*f23
        gun1[:,1] += self.mo_coeff[i,:]*fun1[1]*f13
        gun1[:,2] += self.mo_coeff[i,:]*fun1[2]*f12
         
    rho = numpy.einsum('i,i->', self.mo_occ, gun*gun)
    grad = numpy.zeros(3)
    grad[0] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,0])
    grad[1] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,1])
    grad[2] = numpy.einsum('i,i->', self.mo_occ, gun*gun1[:,2])
    grad = grad + grad
    gradmod = numpy.linalg.norm(grad)

    return rho, grad, gradmod

if __name__ == '__main__':

    from pyscf import lib, dft
    from pypmd.mole import Mole

    log = lib.logger.Logger(sys.stdout, 4)
    log.verbose = 5
    lib.logger.TIMER_LEVEL = 5

    name = 'lif.wfn'
    mol = Mole(name)
    mol.verbose = 4
    mol.kernel()

    a = [0.0, 0.1, 1.0]
    a = numpy.asarray(a)
    t0 = time.clock()
    rho, grad, gradmod = rhograd(a, mol)
    log.timer('simple', t0)
    print a, rho, grad, gradmod

    mol = lib.chkfile.load_mol('lif.chk')
    mf = lib.chkfile.load('lif.chk', 'scf')
    mo_coeff = lib.chkfile.load('lif.chk', 'scf/mo_coeff')
    mo_occ = lib.chkfile.load('lif.chk', 'scf/mo_occ')
    a = numpy.reshape(a, (-1,3))
    t0 = time.clock()
    ao = dft.numint.eval_ao(mol, a, deriv=1)
    rho = dft.numint.eval_rho2(mol, ao, mo_coeff, mo_occ, xctype='GGA')
    log.timer('pyscf', t0)
    gradmod = numpy.linalg.norm(rho[-3:])
    print a, rho, gradmod

