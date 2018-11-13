#!/usr/bin/env python

import numpy, ctypes
from pypmd import lib

libleb = lib.load_library('libsurf')

EPS = 1e-7
LEBEDEV_NGRID = numpy.asarray((
    1   , 6   , 14  , 26  , 38  , 50  , 74  , 86  , 110 , 146 ,
    170 , 194 , 230 , 266 , 302 , 350 , 434 , 590 , 770 , 974 ,
    1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334,
    4802, 5294, 5810))

def lebgrid(npang):

    if npang not in LEBEDEV_NGRID:
        raise ValueError('Lebgrid unsupported angular grid %d' % npang)
    else:
        ct = numpy.zeros((npang), order='F')
        st = numpy.zeros((npang), order='F') 
        cp = numpy.zeros((npang), order='F') 
        sp = numpy.zeros((npang), order='F') 
        aw = numpy.zeros((npang), order='F') 
        agrids = numpy.zeros((npang,5))
        libleb.lebgrid(ct.ctypes.data_as(ctypes.c_void_p), 
                       st.ctypes.data_as(ctypes.c_void_p), 
                       cp.ctypes.data_as(ctypes.c_void_p), 
                       sp.ctypes.data_as(ctypes.c_void_p), 
                       aw.ctypes.data_as(ctypes.c_void_p),  
                       ctypes.c_int(npang)) 
        agrids[:,0] = ct
        agrids[:,1] = st
        agrids[:,2] = cp
        agrids[:,3] = cp
        agrids[:,4] = aw

    return agrids

if __name__ == '__main__':
    npang = 5810
    agrid = lebgrid(npang)
    with open('agrid.txt', 'w') as f2:
        f2.write('# Point 1 2 3 4 weight\n')
        for i in range(npang):
            f2.write('%d   %.6f  %.6f  %.6f  %.6f  %.6f\n' % \
            ((i+1), agrid[i,0], agrid[i,1], agrid[i,2], agrid[i,3], agrid[i,4]))

