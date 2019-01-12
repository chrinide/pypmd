#!/usr/bin/env python

import numpy
import ctypes

import misc

libfapi = misc.load_library('libfapi')

def eval_overlap(self):

    out = numpy.zeros((self.nprims,self.nprims))
    feval = 'overlap'
    drv = getattr(libfapi, feval)
    drv(ctypes.c_int(self.nmo), 
        ctypes.c_int(self.nprims),  
        self.icen.ctypes.data_as(ctypes.c_void_p), 
        self.ityp.ctypes.data_as(ctypes.c_void_p), 
        self.oexp.ctypes.data_as(ctypes.c_void_p), 
        self.ngroup.ctypes.data_as(ctypes.c_void_p), 
        self.nzexp.ctypes.data_as(ctypes.c_void_p), 
        self.nuexp.ctypes.data_as(ctypes.c_void_p), 
        self.mo_coeff.ctypes.data_as(ctypes.c_void_p), 
        self.mo_occ.ctypes.data_as(ctypes.c_void_p), 
        ctypes.c_int(self.natm),  
        self.coords.ctypes.data_as(ctypes.c_void_p), 
        out.ctypes.data_as(ctypes.c_void_p))
        
    return out

def eval_basin(self,coords):

    npoints = coords.shape[0]
    out = numpy.zeros((npoints,3))
    feval = 'eval_basin'
    drv = getattr(libfapi, feval)
    drv(ctypes.c_int(self.nmo), 
        ctypes.c_int(self.nprims),  
        self.icen.ctypes.data_as(ctypes.c_void_p), 
        self.ityp.ctypes.data_as(ctypes.c_void_p), 
        self.oexp.ctypes.data_as(ctypes.c_void_p), 
        self.ngroup.ctypes.data_as(ctypes.c_void_p), 
        self.nzexp.ctypes.data_as(ctypes.c_void_p), 
        self.nuexp.ctypes.data_as(ctypes.c_void_p), 
        self.rcutte.ctypes.data_as(ctypes.c_void_p), 
        self.mo_coeff.ctypes.data_as(ctypes.c_void_p), 
        self.mo_occ.ctypes.data_as(ctypes.c_void_p), 
        ctypes.c_int(self.natm),  
        self.coords.ctypes.data_as(ctypes.c_void_p), 
        ctypes.c_int(npoints),  
        coords.ctypes.data_as(ctypes.c_void_p),
        out.ctypes.data_as(ctypes.c_void_p))
        
    return out.transpose(1,0)

