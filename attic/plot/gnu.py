#!/usr/bin/env python

import h5py
import numpy
from pyscf import dft, lib

def legendre(n):

    x = numpy.zeros(n)
    w = numpy.zeros(n)

    e1 = n*(n+1)
    m = ((n+1)//2)
    for i in range(1,m+1):
        mp1mi = m+1-i
        t = float(4*i-1)*numpy.pi/float(4*n+2)
        x0 = numpy.cos(t)*(1.0-(1.0-1.0/float(n))/float(8*n*n))
        pkm1 = 1.0
        pk = x0
        for k in range(2,n+1):
            pkp1 = 2.0*x0*pk-pkm1-(x0*pk-pkm1)/float(k)
            pkm1 = pk
            pk = pkp1
        d1 = float(n)*(pkm1-x0*pk)
        dpn = d1/(1.0-x0*x0)
        d2pn = (2.0*x0*dpn-e1*pk)/(1.0-x0*x0)
        d3pn = (4.0*x0*d2pn+(2.0-e1)*dpn)/(1.0-x0*x0)
        d4pn = (6.0*x0*d3pn+(6.0-e1)*d2pn)/(1.0-x0*x0)
        u = pk/dpn
        v = d2pn/dpn
        h = -u*(1.0+0.5*u*(v+u*(v*v-d3pn/(3.0*dpn))))
        p = pk+h*(dpn+0.5*h*(d2pn+h/3.0*(d3pn+0.25*h*d4pn)))
        dp = dpn+h*(d2pn+0.5*h*(d3pn+h*d4pn/3.0))
        h = h-p/dp
        xtemp = x0+h
        x[mp1mi-1] = xtemp
        fx = d1-h*e1*(pk+0.5*h*(dpn+h/3.0*(d2pn+0.25*h*(d3pn+0.2*h*d4pn))))
        w[mp1mi-1] = 2.0*(1.0-xtemp*xtemp)/fx/fx

    if ((n%2) == 1):
        x[0] = 0.0
    nmove = ((n+1)//2)
    ncopy = n-nmove
    for i in range(1,nmove+1):
        iback = n+1-i
        x[iback-1] = x[iback-ncopy-1]
        w[iback-1] = w[iback-ncopy-1]
    for i in range(1,n-nmove+1):
        x[i-1] = -x[n-i]
        w[i-1] = w[n-i]

    return x, w

def rquad(nr,r0,rfar,rad,iqudr,mapr):

    rmesh = numpy.zeros(nr)
    dvol = numpy.zeros(nr)
    dvoln = numpy.zeros(nr)
 
    #if (rfar-r0 <= 0.001):
    #    raise RuntimeError('rmax < rmin ??') 
    #if (mapr > 3):
    #    raise RuntimeError('not allowed radial mapping')

    # Determine eta parameter in case of radial mapping
    rfarc = rfar - r0
    if (mapr == 1):
        eta = 2.0*rad/rfarc
    elif (mapr == 2):
        eta = 2.0*numpy.exp(-rfarc/rad)/(1.0-numpy.exp(-rfarc/rad))

    #if (iqudr == 1):
    xr, rwei = legendre(nr)

    # Determine abscissas and volume elements.
    # for finite range (a..b) the transformation is y = (b-a)*x/2+(b+a)/2
    # x = (b-a)*0.5_rp*x+(b+a)*0.5_rp
    # w = w*(b-a)*0.5_rp
    if (mapr == 0):
        for i in range(nr):
            aa = (rfar-r0)/2.0
            bb = (rfar+r0)/2.0
            u = xr[i]
            r = aa*u+bb
            rmesh[i] = r
            dvoln[i] = r*aa
            dvol[i] = dvoln[i]*r
    elif (mapr == 1):
        for i in range(nr):
            u = xr[i]
            den = (1.0-u+eta)
            r = rad*(1.0+u)/den + r0
            rmesh[i] = r
            #if (numpy.abs(den) >= RHOEPS):
            dvoln[i] = rad*(2.0+eta)/den/den*r
            #else
            #dvoln(n) = 0.0_rp
            dvol[i] = dvoln[i]*r
    elif (mapr == 2):
        for i in range(nr):
            u = xr[i]
            den = (1.0-u+eta)
            r = rad*numpy.log((2.0+eta)/den) + r0
            rmesh[i] = r
            dvoln[i] = r*rad/den
            dvol[i] = dvoln[i]*r

    return rmesh, rwei, dvol, dvoln
 

with h5py.File('gamma.chk.h5') as f:
    rmin = f['atom0/rmin'].value
    rmax = f['atom0/rmax'].value
    xyzrho = f['atom0/xyzrho'].value
    xnuc = f['atom0/xnuc'].value
    npang = f['atom0/npang'].value
    ntrial = f['atom0/ntrial'].value
    rsurf = f['atom0/rsurf'].value
    nlimsurf = f['atom0/nlimsurf'].value
    coords = f['atom0/coords'].value

EPS = 1e-6
def inbasin(r,j):

    isin = False
    rs1 = 0.0
    irange = nlimsurf[j]
    irange = int(irange)
    for k in range(irange):
        rs2 = rsurf[j,k]
        if (r >= rs1-EPS and r <= rs2+EPS):
            if (((k+1)%2) == 0):
                isin = False
            else:
                isin = True
            return isin
        rs1 = rs2

    return isin

betafac = 0.5
nrad = 11
iqudr = 1
brad = rmin*betafac
r0 = brad
rfar = rmax
rad = 0.5
mapr = 0
rmesh, rwei, dvol, dvoln = rquad(nrad,r0,rfar,rad,iqudr,mapr)

x = []
y = []
z = []
w = []
xcoor = numpy.zeros(3)
for n in range(nrad):
    r = rmesh[n]
    for j in range(npang):
        inside = True
        inside = inbasin(r,j)
        if (inside == True):
            cost = coords[j,0]
            sintcosp = coords[j,1]*coords[j,2]
            sintsinp = coords[j,1]*coords[j,3]
            xcoor[0] = r*sintcosp
            xcoor[1] = r*sintsinp
            xcoor[2] = r*cost    
            p = xnuc + xcoor
            x.append(p[0])
            y.append(p[1])
            z.append(p[2])
            w.append(coords[j,4]*rwei[n])
            print '%16.8f %16.8f %16.8f %16.8f' % (p[0], p[1], p[2], coords[j,4]*rwei[n])

#x = numpy.array(x)
#y = numpy.array(y)
#z = numpy.array(z)
#w = numpy.array(w)
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#p = ax.scatter(x, y, z, c=w, marker='.')
#plt.colorbar(p)
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')
#plt.show()
