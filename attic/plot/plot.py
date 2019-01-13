#!/usr/bin/env python

import h5py
import numpy, legendre
from pyscf import dft, lib

with h5py.File('gamma.chk.h5') as f:
    print('Surface file is a HDF5 file. data are stored in file/directory structure.')
    print('/', f.keys())
    print('/atom0', f['atom0'].keys())
    rmin = f['atom0/rmin'].value
    rmax = f['atom0/rmax'].value
    xyzrho = f['atom0/xyzrho'].value
    xnuc = f['atom0/xnuc'].value
    npang = f['atom0/npang'].value
    ntrial = f['atom0/ntrial'].value
    rsurf = f['atom0/rsurf'].value
    nlimsurf = f['atom0/nlimsurf'].value
    coords = f['atom0/coords'].value

print 'xnuc', xnuc
print 'xyzrho', xyzrho
print 'npang', npang
print 'ntrial', ntrial
print 'rmin', rmin
print 'rmax', rmax
xdeltain = numpy.zeros(3)
x = numpy.zeros(npang)
y = numpy.zeros(npang)
z = numpy.zeros(npang)
w = numpy.zeros(npang)
points = numpy.zeros((npang,3))
print "ply"
print "format ascii 1.0"
print "comment written by JLC"
print "element vertex 5810"
print "property float x"
print "property float y"
print "property float z"
print "end_header"
for i in range(npang):
    cost = coords[i,0]
    sintcosp = coords[i,1]*coords[i,2]
    sintsinp = coords[i,1]*coords[i,3]
    ract = rsurf[i,0]
    xdeltain[0] = ract*sintcosp
    xdeltain[1] = ract*sintsinp
    xdeltain[2] = ract*cost    
    xpoint = xnuc + xdeltain
    x[i] = xdeltain[0]
    y[i] = xdeltain[1]
    z[i] = xdeltain[2]
    w[i] = coords[i,4]
    points[i] = xpoint
    print '%16.8f %16.8f %16.8f %16.8f' % (xdeltain[0], xdeltain[1], xdeltain[2], coords[i,4])

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = ax.scatter(x, y, z, c=w, marker='o')
plt.colorbar(p)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

