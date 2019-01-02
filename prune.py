#!/usr/bin/env python

NELEC_ERROR_TOL = 0.02
def prune_small_rho_grids(self, grids):
    rho = numint.eval_rho(self, grids.coords)
    n = numpy.dot(rho, grids.weights)
    if abs(n-self.nelectron) < NELEC_ERROR_TOL*n:
        rho *= grids.weights
        idx = abs(rho) > self.small_rho_cutoff/grids.weights.size
        logger.debug(self,'Drop grids %d',
                     grids.weights.size - numpy.count_nonzero(idx))
        grids.coords = numpy.asarray(grids.coords[idx], order='C')
        grids.weights = numpy.asarray(grids.weights[idx], order='C')
    return grids

