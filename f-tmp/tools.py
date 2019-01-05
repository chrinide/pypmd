#!/usr/bin/env python

import h5py
import numpy
import sys

import lib
log = lib.logger.Logger(sys.stdout, 4)

BASE = 0
OUTPUT_COLS = 5
OUTPUT_DIGITS = 5

def dump_tri(stdout, c, ncol=OUTPUT_COLS, digits=OUTPUT_DIGITS, start=BASE):
    ''' Format print for the lower triangular part of an array

    Args:
        stdout : file object
            eg sys.stdout, or stdout = open('/path/to/file') or
            mol.stdout if mol is an object initialized from :class:`gto.Mole`
        c : numpy.ndarray
            coefficients

    Kwargs:
        ncol : int
            Number of columns in the format output (default 5)
        digits : int
            Number of digits of precision for floating point output (default 5)
        start : int
            The number to start to count the index (default 0)

    Examples:

        >>> import sys, numpy
        >>> dm = numpy.eye(3)
        >>> dump_tri(sys.stdout, dm)
                #0        #1        #2   
        0       1.00000
        1       0.00000   1.00000
        2       0.00000   0.00000   1.00000
    '''
    nc = c.shape[1]
    for ic in range(0, nc, ncol):
        dc = c[:,ic:ic+ncol]
        m = dc.shape[1]
        fmt = (' %%%d.%df'%(digits+4,digits))*m + '\n'
        stdout.write(((' '*(digits+3))+'%s\n') % \
                     (' '*(digits)).join(['#%-4d'%i for i in range(start+ic,start+ic+m)]))
        for k, v in enumerate(dc[ic:ic+m]):
            fmt = (' %%%d.%df'%(digits+4,digits))*(k+1) + '\n'
            stdout.write(('%-5d' % (ic+k+start)) + (fmt % tuple(v[:k+1])))
        for k, v in enumerate(dc[ic+m:]):
            stdout.write(('%-5d' % (ic+m+k+start)) + (fmt % tuple(v)))

def print_ply():
    msg = ('Ply format not yet available')
    raise NotImplementedError(msg)

def print_gnu():
    msg = ('Gnuplot format not yet available')
    raise NotImplementedError(msg)

def print_txt(filename, inuc):

    log.info('Surface file is : %s' % (filename))
    idx = 'atom'+str(inuc)
    log.info('Surface info for atom : %d' % inuc)

    with h5py.File(filename) as f:
        i = f[idx+'/inuc'].value
        xnuc = f[idx+'/xnuc'].value
        log.info('Nuclei %d position : %8.5f %8.5f %8.5f', i, *xnuc)
        xyzrho = f[idx+'/xyzrho'].value
        log.info('Nuclei rho %d position : %8.5f %8.5f %8.5f', i, *xyzrho[i])
        npang = f[idx+'/npang'].value
        log.info('Number of angular points : %d', npang)
        ntrial = f[idx+'/ntrial'].value
        log.info('Ntrial : %d', ntrial)
        rsurf = f[idx+'/rsurf'].value
        nlimsurf = f[idx+'/nlimsurf'].value
        coords = f[idx+'/coords'].value
        rmin = f[idx+'/rmin'].value
        rmax = f[idx+'/rmax'].value
        log.info('Rmin and rmax for surface : %8.5f %8.5f', rmin, rmax)
        surf_file = filename+'_'+idx+'.txt'
        with open(surf_file, 'w') as f2:
            for i in range(npang):
                data = str(rsurf[i,:nlimsurf[i]])[1:-1]
                f2.write('%.15f %.15f %.15f %.15f %.15f %s\n' % \
                (coords[i,0],coords[i,1],coords[i,2],coords[i,3],coords[i,4],data))

NPROPS = 3
PROPS = ['density', 'kinetic', 'laplacian']
def print_properties(name, natm):                
    props = numpy.zeros((natm,NPROPS))
    totprops = numpy.zeros((NPROPS))
    with h5py.File(name) as f:
        for i in range(natm):
            idx = 'atom_props'+str(i)
            props[i] = f[idx+'/totprops'].value
    for i in range(natm):
        for j in range(NPROPS):
            log.info('Nuclei %d prop %s value : %8.5f', i, PROPS[j], props[i,j])
            totprops[j] += props[i,j]
    for j in range(NPROPS):
        log.info('Tot prop %s value : %8.5f', PROPS[j], totprops[j])

def read_cube(filename):
    log.info('Cube file is : %s' % (filename))
    with open(filename, 'r') as fin:
        comment1 = fin.readline() #save 1st comment
        comment2 = fin.readline() #save 2nd comment
        norigin = fin.readline().split() # number of atoms and origin
        natoms = int(norigin[0]) #number of atoms
        origin = numpy.array([float(norigin[1]),float(norigin[2]),float(norigin[3])]) #position of origin
        nvoxel = fin.readline().split() #number of voxels
        nx = int(nvoxel[0])
        x = numpy.array([float(nvoxel[1]),float(nvoxel[2]),float(nvoxel[3])])
        nvoxel = fin.readline().split() #
        ny = int(nvoxel[0])
        y = numpy.array([float(nvoxel[1]),float(nvoxel[2]),float(nvoxel[3])])
        nvoxel = fin.readline().split() #
        nz = int(nvoxel[0])
        z = numpy.array([float(nvoxel[1]),float(nvoxel[2]),float(nvoxel[3])])
        atoms = []
        atomsxyz = []
        for atom in range(natoms):
            line= fin.readline().split()
            atoms.append(line[0])
            atomsxyz.append(map(float,[line[2], line[3], line[4]]))
        data = numpy.zeros((nx,ny,nz))
        i = 0
        for s in fin:
            for v in s.split():
               data[i/(ny*nz),(i/nz)%ny,i%nz] = float(v)
               i += 1
        if (i != nx*ny*nz): 
	    raise RuntimeError('wrong number points readed')

    vol = numpy.linalg.det(numpy.array([x,y,z]))
    edensity = numpy.sum(data)
    nelectron = vol*edensity
    log.info('Number of electrons: %.7g' % (nelectron))

#del(BASE,OUTPUT_COLS,OUTPUT_DIGITS,NPROPS,PROPS)

if __name__ == '__main__':
    name = 'h2o.wfn.h5'
    inuc = 0
    print_txt(name,inuc)
    #inuc = 1
    #print_txt(name,inuc)
    #inuc = 2
    #print_txt(name,inuc)

