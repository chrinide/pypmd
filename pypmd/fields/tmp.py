#!/usr/bin/python

def rhograd(point):

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

    return rho, grad, hess, lap, king

