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

import os
import sys
import types
import ctypes
import numpy
import h5py
import tempfile

import param

ASYNC_IO = True
c_double_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int)
c_null_ptr = ctypes.POINTER(ctypes.c_void_p)

if h5py.version.version[:4] == '2.2.':
    sys.stderr.write('h5py-%s is found in your environment. '
                     'h5py-%s has bug in threading mode.\n'
                     'Async-IO is disabled.\n' % ((h5py.version.version,)*2))

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
    libcapi = load_library('libcapi')
    if n is not None:
        libcapi.set_omp_threads.restype = ctypes.c_int
        threads = libcapi.set_omp_threads(ctypes.c_int(int(n)))
        if threads == 0:
            warnings.warn('OpenMP is not available. '
                          'Setting omp_threads to %s has no effects.' % n)
        return threads
    else:
        libcapi.get_omp_threads.restype = ctypes.c_int
        return libcapi.get_omp_threads()

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

class call_in_background(object):
    '''Within this macro, function(s) can be executed asynchronously (the
    given functions are executed in background).

    Attributes:
        sync (bool): Whether to run in synchronized mode.  The default value
            is False (asynchoronized mode).

    Examples:

    >>> with call_in_background(fun) as async_fun:
    ...     async_fun(a, b)  # == fun(a, b)
    ...     do_something_else()

    >>> with call_in_background(fun1, fun2) as (afun1, afun2):
    ...     afun2(a, b)
    ...     do_something_else()
    ...     afun2(a, b)
    ...     do_something_else()
    ...     afun1(a, b)
    ...     do_something_else()
    '''

    def __init__(self, *fns, **kwargs):
        self.fns = fns
        self.handler = None
        self.sync = kwargs.get('sync', not ASYNC_IO)

    if h5py.version.version[:4] == '2.2.': # h5py-2.2.* has bug in threading mode
        # Disable back-ground mode
        def __enter__(self):
            if len(self.fns) == 1:
                return self.fns[0]
            else:
                return self.fns

    else:
        def __enter__(self):
            if self.sync or imp.lock_held():
# Some modules like nosetests, coverage etc
#   python -m unittest test_xxx.py  or  nosetests test_xxx.py
# hang when Python multi-threading was used in the import stage due to (Python
# import lock) bug in the threading module.  See also
# https://github.com/paramiko/paramiko/issues/104
# https://docs.python.org/2/library/threading.html#importing-in-threaded-code
# Disable the asynchoronous mode for safe importing
                def def_async_fn(fn):
                    return fn

            else:
                # Enable back-ground mode
                def def_async_fn(fn):
                    def async_fn(*args, **kwargs):
                        if self.handler is not None:
                            self.handler.join()
                        self.handler = ThreadWithTraceBack(target=fn, args=args,
                                                           kwargs=kwargs)
                        self.handler.start()
                        return self.handler
                    return async_fn

            if len(self.fns) == 1:
                return def_async_fn(self.fns[0])
            else:
                return [def_async_fn(fn) for fn in self.fns]

    def __exit__(self, type, value, traceback):
        if self.handler is not None:
            self.handler.join()


class H5TmpFile(h5py.File):
    '''Create and return an HDF5 temporary file.

    Kwargs:
        filename : str or None
            If a string is given, an HDF5 file of the given filename will be
            created. The temporary file will exist even if the H5TmpFile
            object is released.  If nothing is specified, the HDF5 temporary
            file will be deleted when the H5TmpFile object is released.

    The return object is an h5py.File object. The file will be automatically
    deleted when it is closed or the object is released (unless filename is
    specified).

    Examples:

    >>> from pyscf import lib
    >>> ftmp = lib.H5TmpFile()
    '''
    def __init__(self, filename=None, *args, **kwargs):
        if filename is None:
            tmpfile = tempfile.NamedTemporaryFile(dir=param.TMPDIR)
            filename = tmpfile.name
        h5py.File.__init__(self, filename, *args, **kwargs)
#FIXME: Does GC flush/close the HDF5 file when releasing the resource?
# To make HDF5 file reusable, file has to be closed or flushed
    def __del__(self):
        try:
            self.close()
        except ValueError:  # if close() is called twice
            pass
