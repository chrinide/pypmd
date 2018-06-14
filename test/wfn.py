#!/usr/bin/python

import numpy 

def load_wfn(filename):

    def helper_num(f):
        '''Reads number of orbitals, primitives and atoms'''
        line = f.readline()
        assert line.startswith('GAUSSIAN')
        return [int(i) for i in line.split() if i.isdigit()]

    def helper_coordinates(f):
        '''Reads the coordiantes of the atoms'''
        numbers = numpy.empty(num_atoms, dtype=numpy.int32)
        coordinates = numpy.empty((num_atoms, 3))
        symbols = []
        for atom in range(num_atoms):
            line = f.readline()
            line = line.split()
            symbols.append(line[0])
            numbers[atom] = int(float(line[9]))
            coordinates[atom,:] = [line[4], line[5], line[6]]
        return symbols, numbers, coordinates

    def helper_section(f, start, skip):
        '''Reads CENTRE ASSIGNMENTS, TYPE ASSIGNMENTS, and EXPONENTS sections'''
        section = []
        while len(section) < num_primitives:
            line = f.readline()
            assert line.startswith(start)
            line = line.split()
            section.extend(line[skip:])
        assert len(section) == num_primitives
        return section

    def helper_mo(f):
        '''Reads all MO information'''
        line = f.readline()
        assert line.startswith('MO')
        line = line.split()
        count = line[1]
        occ, energy = line[-5], line[-1]
        coeffs = helper_section(f, ' ', 0)
        coeffs = [i.replace('D', 'E') for i in coeffs]
        return count, occ, energy, coeffs

    with open(filename) as f:
        title = f.readline().strip()
        num_mo, num_primitives, num_atoms = helper_num(f)
        symbols, numbers, coordinates = helper_coordinates(f)
        centers = helper_section(f, 'CENTRE ASSIGNMENTS', 2)
        centers = numpy.array([int(i) for i in centers])
        type_assignment = helper_section(f, 'TYPE ASSIGNMENTS', 2)
        type_assignment = numpy.array([int(i) for i in type_assignment])
        exponents = helper_section(f, 'EXPONENTS', 1)
        exponents = numpy.array([float(i.replace('D', 'E')) for i in exponents])
        mo_count = numpy.empty(num_mo, dtype=numpy.int32)
        mo_occ = numpy.empty(num_mo)
        mo_energy = numpy.empty(num_mo)
        coefficients = numpy.empty([num_primitives, num_mo])
        for mo in range(num_mo):
            mo_count[mo], mo_occ[mo], mo_energy[mo], coefficients[:, mo] = helper_mo(f)
        line = f.readline()

    return title, symbols, numbers, coordinates, centers, type_assignment, exponents, \
        mo_count, mo_occ, mo_energy, coefficients

def read_rdm1(file1, norb):
    f2 = open(file1, 'r')
    rdm1 = numpy.zeros((norb, norb))
    for line in f2.readlines():
        linesp = line.split()
        ind1 = int(linesp[0])
        ind2 = int(linesp[1])
        ind1 = ind1 - 1
        ind2 = ind2 - 1
        rdm1[ind1, ind2] = float(linesp[2])

    f2.close()

    return rdm1

if __name__ == '__main__':

    import data # Almacenas los tipos de primitivas
    name = 'n2_rhf.wfn'

    title, symbols, numbers, coords, icen, ityp, oexp, \
    mo_count, mo_occ, mo_energy, mo_coeff = load_wfn(name)

    nprims, nmo = mo_coeff.shape

    #rdm1 = read_rdm1('rdm1.dat', nmo)
    #print "One particle density matrix is"
    #print rdm1
    #
    #natocc, natorb = numpy.linalg.eigh(-rdm1)
    #for i, k in enumerate(numpy.argmax(abs(natorb), axis=0)):
    #    if natorb[k,i] < 0:
    #        natorb[:,i] *= -1
    #natorb = numpy.dot(mo_coeff[:,:], natorb)
    #natocc = -natocc

    def density(point):
        fun = numpy.zeros((3))
        fun1 = numpy.zeros((3))
        fun2 = numpy.zeros((3))
        gun = numpy.zeros((nmo))
        gun1 = numpy.zeros((nmo,3))
        gun2 = numpy.zeros((nmo,6))
        xcoor = numpy.zeros((3))
        for i in range(nprims):
            ic = icen[i]-1
            itip = ityp[i]-1
            it = data.nlm[itip,:]
            ori = -oexp[i]
            dp2 = ori+ori
            xcoor = point - coords[ic,:]
            dis2 = numpy.einsum('i,i->', xcoor, xcoor)
            aexp = numpy.exp(ori*dis2)
            for j in range(3):
                n = it[j]
                x = xcoor[j]
                if (n == 0):
                    dp2x = dp2*x
                    dp2x2 = dp2*x*x
                    fun2[j] = dp2*(1.0+dp2x2)
                    fun1[j] = dp2x
                    fun[j] = 1.0
                elif (n == 1):
                    x2 = x*x
                    dp2x2 = dp2*x2
                    fun2[j] = dp2*x*(3.0+dp2x2)
                    fun1[j] = 1.0+dp2x2
                    fun[j] = x
                elif (n == 2):
                    x2 = x*x
                    dp2x2 = dp2*x2
                    fun2[j] = 2.0+dp2x2*(5.0+dp2x2)
                    fun1[j] = x*(2.0+dp2x2)
                    fun[j] = x2
                elif (n == 3):
                    x2 = x*x
                    dp2x2 = dp2*x2
                    fun2[j] = x*(6.0+dp2x2*(7.0+dp2x2))
                    fun1[j] = x2*(3.0+dp2x2)
                    fun[j] = x*x2
                elif (n == 4):
                    x2 = x*x
                    dp2x2 = dp2*x2
                    fun2[j] = x2*(x2*(dp2*(dp2x2+9.0))+12.0)
                    fun1[j] = x2*x*(4.0+dp2x2)
                    fun[j] = x2*x2
                elif (n == 5):
                    x2 = x*x
                    dp2x2 = dp2*x2
                    fun2[j] = x2*x*(x2*(dp2*(dp2x2+11.0))+20.0) 
                    fun1[j] = x2*x2*(5.0+dp2x2)
                    fun[j] = x2*x2*x

            f23 = fun[1]*fun[2]*aexp
            f13 = fun[0]*fun[2]*aexp
            f12 = fun[0]*fun[1]*aexp
            g23 = fun1[1]*fun[2]*aexp
            g32 = fun1[2]*fun[1]*aexp
            g21 = fun1[1]*fun[0]*aexp
            for j in range(nmo):
                cfj = mo_coeff[i,j]
                gun[j] += cfj*fun[0]*f23
                gun1[j,0] += cfj*(fun1[0]*f23)
                gun1[j,1] += cfj*(fun1[1]*f13)
                gun1[j,2] += cfj*(fun1[2]*f12)
                gun2[j,0] += cfj*(fun2[0]*f23)
                gun2[j,2] += cfj*(fun2[1]*f13)
                gun2[j,5] += cfj*(fun2[2]*f12)
                gun2[j,1] += cfj*(fun1[0]*g23)
                gun2[j,3] += cfj*(fun1[0]*g32)
                gun2[j,4] += cfj*(fun1[2]*g21)
             
        rho = numpy.einsum('i,i->', mo_occ, gun*gun)
        grad = numpy.zeros(3)
        grad[0] = numpy.einsum('i,i->', mo_occ, gun*gun1[:,0])
        grad[1] = numpy.einsum('i,i->', mo_occ, gun*gun1[:,1])
        grad[2] = numpy.einsum('i,i->', mo_occ, gun*gun1[:,2])
        grad = grad + grad
        gradmod = numpy.linalg.norm(grad)
        king = numpy.einsum('i,i->', mo_occ, gun1[:,0]*gun1[:,0])
        king += numpy.einsum('i,i->', mo_occ, gun1[:,1]*gun1[:,1])
        king += numpy.einsum('i,i->', mo_occ, gun1[:,2]*gun1[:,2])
        king = king*0.5
        hess = numpy.zeros((3,3))
        hess[0,0] = numpy.einsum('i,i->', mo_occ, gun2[:,0]*gun)
        hess[0,0] += numpy.einsum('i,i->', mo_occ, gun1[:,0]*gun1[:,0]) 
        hess[0,0] = hess[0,0] + hess[0,0]
        hess[1,1] = numpy.einsum('i,i->', mo_occ, gun2[:,2]*gun)
        hess[1,1] += numpy.einsum('i,i->', mo_occ,gun1[:,1]*gun1[:,1]) 
        hess[1,1] = hess[1,1] + hess[1,1]
        hess[2,2] = numpy.einsum('i,i->', mo_occ, gun2[:,5]*gun)
        hess[2,2] += numpy.einsum('i,i->', mo_occ, gun1[:,2]*gun1[:,2]) 
        hess[2,2] = hess[2,2] + hess[2,2]
        lap = hess[0,0] + hess[1,1] + hess[2,2]

        hess[0,1] = numpy.einsum('i,i->', mo_occ, gun2[:,1]*gun)
        hess[0,1] += numpy.einsum('i,i->', mo_occ, gun1[:,0]*gun1[:,1])
        hess[0,1] = hess[0,1] + hess[0,1]
        hess[1,0] = hess[0,1]

        hess[0,2] = numpy.einsum('i,i->', mo_occ, gun2[:,3]*gun)
        hess[0,2] += numpy.einsum('i,i->', mo_occ, gun1[:,0]*gun1[:,2])
        hess[0,2] = hess[0,2] + hess[0,2]
        hess[2,0] = hess[0,2]

        hess[1,2] = numpy.einsum('i,i->', mo_occ, gun2[:,4]*gun)
        hess[1,2] += numpy.einsum('i,i->', mo_occ, gun1[:,1]*gun1[:,2])
        hess[1,2] = hess[1,2] + hess[1,2]
        hess[2,1] = hess[1,2]

        return rho, grad, gradmod, hess, lap, king

    point = numpy.zeros(3)
    rho,grad,gradmod,hess,lap,king = density(point)
    print rho, grad, gradmod
    print hess, lap
    print king
