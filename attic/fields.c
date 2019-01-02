#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "mole.h"
#include "basis.h"
#include "fields.h"

void driver(const int nmo, 
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
            const double *coords){

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

  double point[3],grad[3];
	double rho, gradmod;
  point[0] = 0.0;
  point[1] = 0.0;
  point[2] = 0.3;
  point_rho_grad(point, &rho, grad, &gradmod);
  printf("The value of rho is %f\n", rho);
  printf("The value of grad is %f %f %f\n", grad[0], grad[1], grad[2]);
  printf("The value of gradmod is %f\n", gradmod);

	// Pointers for shell version 
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

  point_rho_grad_shell(point, &rho, grad, &gradmod);
  printf("The value of shell rho is %f\n", rho);
  printf("The value of shell grad is %f %f %f\n", grad[0], grad[1], grad[2]);
  printf("The value of shell gradmod is %f\n", gradmod);

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

}

inline void point_rho_grad_shell(const double *p, double *rho, double *grad, double *gradmod)
{
  double gun[nmo_], gun1[nmo_][3];
  double xcoor[3], fun[3], fun1[3];
  int ii, ic, m, jj, j, it[3];

  *rho = 0.0;
  *gradmod = 0.0;
  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = 0.0;
  for (ii=0; ii<nmo_; ii++) {
    gun[ii] = 0.0;
    gun1[ii][0] = 0.0;
    gun1[ii][1] = 0.0;
    gun1[ii][2] = 0.0;
  }
    
	// TODO: avoid natm loop and loop only over total shells
  for (ic=0; ic<natm_; ic++){
    xcoor[0] = p[0] - coords_[ic*3+0];
    xcoor[1] = p[1] - coords_[ic*3+1];
    xcoor[2] = p[2] - coords_[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
  	for (m=0; m<ngroup_[ic]; m++){
    	if (dis2 >= rcutte_[m*natm_+ic]){continue;}
			int idx1 = m*(ngtoh*natm_) + 0*natm_ + ic;
      int k = nuexp_[idx1];
      double ori = -oexp_[k-1];
      double dp2 = ori + ori;
      double aexp = exp(ori*dis2);
      for (jj=0; jj<nzexp_[m*natm_+ic]; jj++){
				int idx2 = m*(ngtoh*natm_) + jj*natm_ + ic;
      	int i = nuexp_[idx2];
        int itip = ityp_[i-1]-1;
    		it[0] = nlm[itip][0];
		    it[1] = nlm[itip][1];
		    it[2] = nlm[itip][2];
    		for (j=0; j<3; j++) {
		      int n = it[j];
		      double x = xcoor[j];
		      if (n==0) {
		        double dp2x = dp2*x;
		        fun1[j] = dp2x;
		        fun[j] = 1.0;
		      }
		      else if (n==1) {
		        double x2 = x*x;
		        double dp2x2 = dp2*x2;
		        fun1[j] = 1.0+dp2x2;
		        fun[j] = x;
		      }
		      else if (n==2) {
		        double x2 = x*x;
		        double dp2x2 = dp2*x2;
		        fun1[j] = x*(2.0+dp2x2);
		        fun[j] = x2;
		      }
		      else if (n==3) {
		        double x2 = x*x;
		        double dp2x2 = dp2*x2;
		        fun1[j] = x2*(3.0+dp2x2);
		        fun[j] = x*x2;
		      }
		      else if (n==4) {
		        double x2 = x*x;
		        double dp2x2 = dp2*x2;
		        fun1[j] = x2*x*(4.0+dp2x2);
		        fun[j] = x2*x2;
		      }
		      else if (n==5) {
		        double x2 = x*x;
		        double dp2x2 = dp2*x2;
		        fun1[j] = x2*x2*(5.0+dp2x2);
		        fun[j] = x2*x2*x;
		      }
		    }
        double f12 = fun[0]*fun[1]*aexp;
        double f123 = f12*fun[2];
        double fa = fun1[0]*fun[1]*fun[2]*aexp;
        double fb = fun1[1]*fun[0]*fun[2]*aexp;
        double fc = fun1[2]*f12;
    		for (j=0; j<nmo_; j++) {
      		double cfj = mo_coeff_[(i-1)*nmo_+j];
          gun[j] += cfj*f123;
          gun1[j][0] += cfj*fa;
          gun1[j][1] += cfj*fb;
          gun1[j][2] += cfj*fc;
				}
      }
    }
  }

  // Run again over orbitals
  for (ii=0; ii<nmo_; ii++) {
    *rho += mo_occ_[ii]*gun[ii]*gun[ii];
    grad[0] += mo_occ_[ii]*gun[ii]*gun1[ii][0];
    grad[1] += mo_occ_[ii]*gun[ii]*gun1[ii][1];
    grad[2] += mo_occ_[ii]*gun[ii]*gun1[ii][2];
  }
  grad[0] += grad[0];
  grad[1] += grad[1];
  grad[2] += grad[2];
  *gradmod = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);

}

// Slow versions no shells
inline void point_rho_grad(const double *p, double *rho, double *grad, double *gradmod)
{

  double gun[nmo_], gun1[nmo_][3];
  double xcoor[3], fun[3], fun1[3];
  int i, j, it[3];

  *rho = 0.0;
  *gradmod = 0.0;
  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = 0.0;
  for (i=0; i<nmo_; i++) {
    gun[i] = 0.0;
    gun1[i][0] = 0.0;
    gun1[i][1] = 0.0;
    gun1[i][2] = 0.0;
  }

  // Simpler version over primitive, no screaning, no cuttof !!
  for (i=0; i<nprims_; i++) {
    int ic = icen_[i]-1;
    int itip = ityp_[i]-1;
    it[0] = nlm[itip][0];
    it[1] = nlm[itip][1];
    it[2] = nlm[itip][2];
    double ori = -oexp_[i];
    double dp2 = ori+ori;
    xcoor[0] = p[0] - coords_[ic*3+0];
    xcoor[1] = p[1] - coords_[ic*3+1];
    xcoor[2] = p[2] - coords_[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
    double aexp = exp(ori*dis2);
    for (j=0; j<3; j++) {
      int n = it[j];
      double x = xcoor[j];
      if (n==0) {
        double dp2x = dp2*x;
        fun1[j] = dp2x;
        fun[j] = 1.0;
      }
      else if (n==1) {
        double x2 = x*x;
        double dp2x2 = dp2*x2;
        fun1[j] = 1.0+dp2x2;
        fun[j] = x;
      }
      else if (n==2) {
        double x2 = x*x;
        double dp2x2 = dp2*x2;
        fun1[j] = x*(2.0+dp2x2);
        fun[j] = x2;
      }
      else if (n==3) {
        double x2 = x*x;
        double dp2x2 = dp2*x2;
        fun1[j] = x2*(3.0+dp2x2);
        fun[j] = x*x2;
      }
      else if (n==4) {
        double x2 = x*x;
        double dp2x2 = dp2*x2;
        fun1[j] = x2*x*(4.0+dp2x2);
        fun[j] = x2*x2;
      }
      else if (n==5) {
        double x2 = x*x;
        double dp2x2 = dp2*x2;
        fun1[j] = x2*x2*(5.0+dp2x2);
        fun[j] = x2*x2*x;
      }
    }
    double f12 = fun[0]*fun[1]*aexp;
    double f123 = fun[2]*f12;
    double fa = fun1[0]*fun[1]*fun[2]*aexp;
    double fb = fun1[1]*fun[0]*fun[2]*aexp;
    double fc = fun1[2]*f12;
    for (j=0; j<nmo_; j++) {
      double cfj = mo_coeff_[i*nmo_+j];
      gun[j] += cfj*f123;
      gun1[j][0] += cfj*fa;
      gun1[j][1] += cfj*fb;
      gun1[j][2] += cfj*fc;
    }
  }

  // Run again over orbitals
  for (i=0; i<nmo_; i++) {
    *rho += mo_occ_[i]*gun[i]*gun[i];
    grad[0] += mo_occ_[i]*gun[i]*gun1[i][0];
    grad[1] += mo_occ_[i]*gun[i]*gun1[i][1];
    grad[2] += mo_occ_[i]*gun[i]*gun1[i][2];
  }
  grad[0] += grad[0];
  grad[1] += grad[1];
  grad[2] += grad[2];
  *gradmod = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);

}

