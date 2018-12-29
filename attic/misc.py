#!/usr/bin/env python
# Copyright 2014-2018 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
Some hacky functions
'''

import os, sys
import types
import ctypes
import numpy

c_double_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int)
c_null_ptr = ctypes.POINTER(ctypes.c_void_p)

def load_library(libname):
# numpy 1.6 has bug in ctypeslib.load_library, see numpy/distutils/misc_util.py
    if '1.6' in numpy.__version__:
        if (sys.platform.startswith('linux') or
            sys.platform.startswith('gnukfreebsd')):
            so_ext = '.so'
        elif sys.platform.startswith('darwin'):
            so_ext = '.dylib'
        elif sys.platform.startswith('win'):
            so_ext = '.dll'
        else:
            raise OSError('Unknown platform')
        libname_so = libname + so_ext
        return ctypes.CDLL(os.path.join(os.path.dirname(__file__), libname_so))
    else:
        _loaderpath = os.path.dirname(__file__)
        return numpy.ctypeslib.load_library(libname, _loaderpath)

def num_threads(n=None):
    '''Set the number of OMP threads.  If argument is not specified, the
    function will return the total number of available OMP threads.

    Examples:

    >>> from pyscf import lib
    >>> print(lib.num_threads())
    8
    >>> lib.num_threads(4)
    4
    >>> print(lib.num_threads())
    4
    '''
    from pyscf.lib.numpy_helper import _np_helper
    if n is not None:
        _np_helper.set_omp_threads.restype = ctypes.c_int
        threads = _np_helper.set_omp_threads(ctypes.c_int(int(n)))
        if threads == 0:
            warnings.warn('OpenMP is not available. '
                          'Setting omp_threads to %s has no effects.' % n)
        return threads
    else:
        _np_helper.get_omp_threads.restype = ctypes.c_int
        return _np_helper.get_omp_threads()

class with_omp_threads(object):
    '''Using this macro to create a temporary context in which the number of
    OpenMP threads are set to the required value. When the program exits the
    context, the number OpenMP threads will be restored.

    Args:
        nthreads : int

    Examples:

    >>> from pyscf import lib
    >>> print(lib.num_threads())
    8
    >>> with lib.with_omp_threads(2):
    ...     print(lib.num_threads())
    2
    >>> print(lib.num_threads())
    8
    '''
    def __init__(self, nthreads=None):
        self.nthreads = nthreads
        self.sys_threads = None
    def __enter__(self):
        if self.nthreads is not None and self.nthreads >= 1:
            self.sys_threads = num_threads()
            num_threads(self.nthreads)
        return self
    def __exit__(self, type, value, traceback):
        if self.sys_threads is not None:
            num_threads(self.sys_threads)


