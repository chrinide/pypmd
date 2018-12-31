#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "mole.h"
#include "basis.h"
#include "vv10_gpu.h"

extern "C" {
void driver_vv10_gpu(const int nmo, 
                 const int nprims,
                 const int *icen, 
                 const int *ityp,
                 const double *oexp, 
                 const int *ngroup,
            		 const int *nzexp,
            		 const int *nuexp,
            		 const double *rcutte,
            		 const double *mo_coeff,
            		 const double *mo_occ, 
            		 const int natm,
            		 const double *coords,
								 const int npoints,
								 const double *points,
								 double *rho,
								 double *gnorm){

  int i, j, k;

  init_nlm();

  // Basis info
  natm_ = natm;
  nprims_ = nprims;
  nmo_ = nmo;
  mo_coeff_ = (double *) malloc(sizeof(double)*nmo_*nprims_);
  assert(mo_coeff_ != NULL);
  mo_occ_ = (double *) malloc(sizeof(double)*nmo_);
  assert(mo_occ_ != NULL);
  icen_ = (int *) malloc(sizeof(int)*nprims_);
  assert(icen_ != NULL);
  ityp_ = (int *) malloc(sizeof(int)*nprims_);
  assert(ityp_ != NULL);
  oexp_ = (double *) malloc(sizeof(double)*nprims_);
  assert(oexp_ != NULL);
  for (i=0; i<nprims_; i++) {
    icen_[i] = icen[i];
    ityp_[i] = ityp[i];
    oexp_[i] = oexp[i];
  }
  for (i=0; i<nmo_; i++) {                      // Orbital
    mo_occ_[i] = mo_occ[i];
    for (j=0; j<nprims_; j++) {
      mo_coeff_[j*nmo_+i] = mo_coeff[j*nmo_+i];
    }
  }
  coords_ = (double *) malloc(sizeof(double)*3*natm_);
  assert(coords_ != NULL);
  for (i=0; i<natm_; i++) {
    coords_[i*3+0] = coords[i*3+0];
    coords_[i*3+1] = coords[i*3+1];
    coords_[i*3+2] = coords[i*3+2];
  }

  ngroup_ = (int *) malloc(sizeof(int)*natm_);
  assert(ngroup_ != NULL);
  for (i=0; i<natm_; i++) {
		ngroup_[i] = ngroup[i];
	}
  nzexp_ = (int *) malloc(sizeof(int)*natm_*mgrp);
  assert(nzexp_ != NULL);
  rcutte_ = (double *) malloc(sizeof(double)*natm_*mgrp);
  assert(rcutte_ != NULL);
  for (i=0; i<natm_; i++) {
    for (j=0; j<ngroup[i]; j++) {
      rcutte_[j*natm_+i] = rcutte[j*natm_+i];
      nzexp_[j*natm_+i] = nzexp[j*natm_+i];
		}
	}
  nuexp_ = (int *) malloc(sizeof(int)*natm_*mgrp*ngtoh);
  assert(nuexp_ != NULL);
  for (i=0; i<natm_; i++) {
    for (j=0; j<ngroup_[i]; j++) {
      for (k=0; k<nzexp_[j*natm_+i]; k++) {
				int idx = j*(ngtoh*natm_) + k*natm_ + i;
        nuexp_[idx] = nuexp[idx];
      }
    }
  }
    
  const int blockSize = 1024;
  const int gridSize = (npoints + blockSize - 1) / blockSize;

  cudaMemcpyToSymbol(dnprims_, &nprims_, sizeof(int));
  cudaMemcpyToSymbol(dnmo_, &nmo_, sizeof(int));
  cudaMemcpyToSymbol(dnatm_, &natm_, sizeof(int));
  cudaMemcpyToSymbol(dnlm_, &nlm, 56*3*sizeof(int));

  free(rcutte_);
	free(ngroup_);
	free(nzexp_);
  free(nuexp_);
  free(mo_coeff_);
  free(mo_occ_);
  free(oexp_);
  free(coords_);
  free(icen_);
  free(ityp_);

}}
