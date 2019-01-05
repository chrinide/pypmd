#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "csurf.h"

void csurf_driver(const int nmo, 
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
                  const int npang,
                  const int inuc,
                  const double *xyzrho,
                  const double *ct,
                  const double *st,
                  const double *cp,
                  const double *sp,
                  const int backend,
                  const int ntrial,
                  const double epsiscp,
                  const double epsroot,
                  const double rmaxsurf,
                  const double epsilon,
                  const double *rpru,
                  const double step,
                  const int mstep,
                  int *nlimsurf, double *rlimsurf){

  int i, j, k;

  // Sizes
  ncent_ = natm;
  nprims_ = nprims;
  nmo_ = nmo;
  inuc_ = inuc;
  epsiscp_ = epsiscp;
  ntrial_ = ntrial;      
  npang_ = npang;
  epsroot_ = epsroot;
  rmaxsurf_ = rmaxsurf;
  backend_ = backend;
  epsilon_ = epsilon;
  step_ = step;
  mstep_ = mstep;
  steeper_ = backend;

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
  coords_ = (double *) malloc(sizeof(double)*3*ncent_);
  assert(coords_ != NULL);
  for (i=0; i<ncent_; i++) {
    coords_[i*3+0] = coords[i*3+0];
    coords_[i*3+1] = coords[i*3+1];
    coords_[i*3+2] = coords[i*3+2];
  }
  xyzrho_ = (double *) malloc(sizeof(double)*3*ncent_);
  assert(xyzrho_ != NULL);
  for (i=0; i<ncent_; i++) {
    xyzrho_[i*3+0] = xyzrho[i*3+0];
    xyzrho_[i*3+1] = xyzrho[i*3+1];
    xyzrho_[i*3+2] = xyzrho[i*3+2];
  }
  xnuc_[0] = xyzrho[inuc_*3+0];
  xnuc_[1] = xyzrho[inuc_*3+1]; 
  xnuc_[2] = xyzrho[inuc_*3+2];

	// Pointers for shell version 
  ngroup_ = (int *) malloc(sizeof(int)*ncent_);
  assert(ngroup_ != NULL);
  for (i=0; i<ncent_; i++) {
		ngroup_[i] = ngroup[i];
	}
  nzexp_ = (int *) malloc(sizeof(int)*ncent_*mgrp);
  assert(nzexp_ != NULL);
  rcutte_ = (double *) malloc(sizeof(double)*ncent_*mgrp);
  assert(rcutte_ != NULL);
  for (i=0; i<ncent_; i++) {
    for (j=0; j<ngroup[i]; j++) {
      rcutte_[j*ncent_+i] = rcutte[j*ncent_+i];
      nzexp_[j*ncent_+i] = nzexp[j*ncent_+i];
		}
	}
  nuexp_ = (int *) malloc(sizeof(int)*ncent_*mgrp*ngtoh);
  assert(nuexp_ != NULL);
  for (i=0; i<ncent_; i++) {
    for (j=0; j<ngroup_[i]; j++) {
      for (k=0; k<nzexp_[j*ncent_+i]; k++) {
				int idx = j*(ngtoh*ncent_) + k*ncent_ + i;
        nuexp_[idx] = nuexp[idx];
      }
    }
  }

	ct_ = (double *) malloc(sizeof(double)*npang_);
  assert(ct_ != NULL);
  st_ = (double *) malloc(sizeof(double)*npang_);
  assert(st_ != NULL);
	cp_ = (double *) malloc(sizeof(double)*npang_);
  assert(cp_ != NULL);
	sp_ = (double *) malloc(sizeof(double)*npang_);
  assert(sp_ != NULL);
  for (i=0; i<npang_; i++){
    ct_[i] = ct[i];
    st_[i] = st[i];
    cp_[i] = cp[i];
    sp_[i] = sp[i];
  	//printf("ct st cp sp %f %f %f %f\n", ct_[i], st_[i], cp_[i], sp_[i]);
  }

	rpru_ = (double *) malloc(sizeof(double)*ntrial_);
  assert(rpru_ != NULL);
  for (i=0; i<ntrial_; i++){
    rpru_[i] = rpru[i];
  }
	rsurf_ = (double *) malloc(sizeof(double)*npang_*ntrial_);
  assert(rsurf_ != NULL);
	nlimsurf_ = (int *) malloc(sizeof(int)*npang_);
  assert(nlimsurf_ != NULL);

  init_nlm();
  surface();
	for (i=0; i<npang_; i++){
    nlimsurf[i] = nlimsurf_[i];
	  for (j=0; j<ntrial_; j++){
      rlimsurf[i*ntrial_+j] = rsurf_[i*ntrial_+j];
    }
  }

  free(rsurf_);
  free(nlimsurf_);
  free(rpru_);
  free(rcutte_);
	free(ngroup_);
	free(nzexp_);
  free(nuexp_);
  free(mo_coeff_);
  free(mo_occ_);
  free(oexp_);
  free(coords_);
  free(xyzrho_);
  free(icen_);
  free(ityp_);
  free(ct_);
  free(st_);
  free(cp_);
  free(sp_);

}

/*
int *mat = (int *)malloc(rows * cols * sizeof(int));
Then, you simulate the matrix using
int offset = i * cols + j;
// now mat[offset] corresponds to m(i, j)
for row-major ordering and
int offset = i + rows * j;
// not mat[offset] corresponds to m(i, j)
*/

void surface(){

  int i, j;
  double xsurf[ntrial_][3]; 
  int isurf[ntrial_][2];
  int nintersec = 0;
  double xpoint[3];

  if (ncent_ == 1){  
	  for (i=0; i<npang_; i++){
      nlimsurf_[i] = 1;
      int offset = i*ntrial_;
      rsurf_[offset] = rmaxsurf_;
      return;
    }
  }

#pragma omp parallel default(none)  \
    private(i,nintersec,j,xpoint,xsurf,isurf) \
    shared(npang_,ct_,st_,cp_,sp_,xnuc_,inuc_,rpru_,step_,\
    ntrial_,epsroot_,rsurf_,nlimsurf_,rmaxsurf_,epsilon_)
{
#pragma omp for schedule(dynamic) nowait
	for (i=0; i<npang_; i++){
    nintersec = 0;
    double cost = ct_[i];
    double sintcosp = st_[i]*cp_[i];
    double sintsinp = st_[i]*sp_[i];
    int ia = inuc_, ib;
    double ra = 0.0, rb;
	  for (j=0; j<ntrial_; j++){
      double ract = rpru_[j];
      xpoint[0] = xnuc_[0] + ract*sintcosp; 
      xpoint[1] = xnuc_[1] + ract*sintsinp; 
      xpoint[2] = xnuc_[2] + ract*cost;     
      //TODO: Better Check for error
      int ier = odeint(xpoint, step_, epsilon_);
      if (ier == 1) {
        cerror("too short steep on odeint");
      } else if (ier == 2) {
        cerror("too many steeps on odeint");
			} else if (ier == 4) {
        cerror("nna on odeint");
			}
      bool good = checkcp(xpoint, &ib);
      rb = ract;
      if (ib != ia && (ia == inuc_ || ib == inuc_)){
        if (ia != inuc_ || ib != -1){
          nintersec += 1;
          xsurf[nintersec-1][0] = ra;
          xsurf[nintersec-1][1] = rb;
          isurf[nintersec-1][0] = ia;
          isurf[nintersec-1][1] = ib;
        }                             
      }
      ia = ib;
      ra = rb;
    }
	  for (j=0; j<nintersec; j++){
      ia = isurf[j][0];
      ib = isurf[j][1];
      ra = xsurf[j][0];
      rb = xsurf[j][1];
      double xin[3], xfin[3];
      xin[0] = xnuc_[0] + ra*sintcosp;
      xin[1] = xnuc_[1] + ra*sintsinp;
      xin[2] = xnuc_[2] + ra*cost;
      xfin[0] = xnuc_[0] + rb*sintcosp;
      xfin[1] = xnuc_[1] + rb*sintsinp;
      xfin[2] = xnuc_[2] + rb*cost;
      while (fabs(ra-rb) > epsroot_){
        double xmed[3];
        xmed[0] = 0.5*(xfin[0] + xin[0]);   
        xmed[1] = 0.5*(xfin[1] + xin[1]);   
        xmed[2] = 0.5*(xfin[2] + xin[2]);   
        double rm = 0.5*(ra + rb);
        xpoint[0] = xmed[0];
        xpoint[1] = xmed[1];
        xpoint[2] = xmed[2];
        int im;
        int ier = odeint(xpoint, step_, epsilon_);
      	if (ier == 1) {
        	cerror("too short steep on odeint");
	      } else if (ier == 2) {
        	cerror("too many steeps on odeint");
				} else if (ier == 4) {
        	cerror("nna on odeint");
				}
        bool good = checkcp(xpoint, &im);
        if (im == ia){
          xin[0] = xmed[0];
          xin[1] = xmed[1];
          xin[2] = xmed[2];
          ra = rm;
        }
        else if (im == ib){
          xfin[0] = xmed[0];
          xfin[1] = xmed[1];
          xfin[2] = xmed[2];
          rb = rm;
        }
        else{
          if (ia == inuc_){
            xfin[0] = xmed[0];
            xfin[1] = xmed[1];
            xfin[2] = xmed[2];
            rb = rm;
          }
          else{
            xin[0] = xmed[0];
            xin[1] = xmed[1];
            xin[2] = xmed[2];
            ra = rm;
          }
        }
      }
      xsurf[j][2] = 0.5*(ra + rb);
    }
    // organize pairs
    nlimsurf_[i] = nintersec; 
	  for (j=0; j<nintersec; j++){
      int offset = i*ntrial_+j;
      rsurf_[offset] = xsurf[j][2];
    }
    if (nintersec%2 == 0){
      nintersec += 1;
      nlimsurf_[i] = nintersec;
      int offset = i*ntrial_+(nintersec-1);
      rsurf_[offset] = rmaxsurf_;
    }
    //printf("#* %d %d %.6f %.6f %.6f %.6f ",i,nlimsurf_[i],ct_[i],st_[i],cp_[i],sp_[i]);
	  //for (j=0; j<nlimsurf_[i]; j++){
    // printf(" %.6f ",rsurf_[i*ntrial_+j]);
    //}
    //printf("\n");
  }
}

}

//ier = 0 (correct), 1 (short step), 2 (too many iterations), 
//      3 (infty), 4 (nna), 5(undef)
int odeint(double *ystart, double h1, double eps){

	double rho, gradmod, hnext;
  double grad[3] = {0.0};
	int ier = 0, i, nuc;
  double dydx[3], y[3], yscal[3];

	rho_grad(ystart, &rho, grad, &gradmod);
  if (rho <= RHOEPS && gradmod <= GRADEPS){
    ier = 3;
    return ier;
  }
  
  double hmin = 0.0;
  double x1 = 0.0;
  double x2 = 1e40;
  double x = x1;
  double h = h1;
  y[0] = ystart[0];
  y[1] = ystart[1];
  y[2] = ystart[2];

	for (i=0; i<mstep_; i++){
	  rho_grad(y, &rho, dydx, &gradmod);
		yscal[0] = fmax(fabs(y[0]) + fabs(dydx[0]*h) + TINY, eps);
		yscal[1] = fmax(fabs(y[1]) + fabs(dydx[1]*h) + TINY, eps);
		yscal[2] = fmax(fabs(y[2]) + fabs(dydx[2]*h) + TINY, eps);
		if ((x+h-x2)*(x+h-x1) > 0.0) h = x2 - x;
		rkqs(y, dydx, &x, h, eps, yscal, &hnext);
		if ((x-x2)*(x2-x1) >= 0.0 || checkcp(y, &nuc)){
			ystart[0] = y[0];
			ystart[1] = y[1];
			ystart[2] = y[2];
      ier = 0;
			return ier;
		}
		if (fabs(hnext) <= hmin) cerror("Step size too small in odeint");
		if (i == (mstep_-1)) cerror("Reached max steps in odeint");
		h = hnext;
  }
    
	// Test if the point is far from RMAXSURF from current atom. 
  double a1 = y[0] - xnuc_[0];
  double a2 = y[1] - xnuc_[1];
  double a3 = y[2] - xnuc_[2];
  if ((a1*a1+a2*a2+a3*a3)>=5.0*5.0){
    ier = 3;
  } else { 
  	printf("NNA at %f %f %f\n", y[0], y[1], y[2]);
	  cerror("NNA found in odeint"); 
  }

  return ier;
}

void rkqs(double *y, double *dydx, double *x, 
          double htry, double eps,
	        double *yscal, double *hnext){
	
  int i;
	double yerr[3], ytemp[3], errmax, xnew, htemp;

	double h = htry;
  *hnext = 0.0;
  errmax = 0.0;

	for (;;){
		steeper_rkck(y, dydx, h, ytemp, yerr);
		errmax = 0.0;
		for (i=0; i<3; i++) {
		  errmax = fmax(errmax, fabs(yerr[i]/yscal[i]));
    }
		errmax /= eps;
		if (errmax > 1.0) {
			htemp = SAFETY*h*pow(errmax, PSHRNK);
      h = fmin(fmax(fabs(htemp),0.1*fabs(h)),h);
			xnew = *x + h;
			if (xnew == *x) {
        cerror("stepsize underflow in rkqs");
      }
			continue;
		}
		else {
			if (errmax > ERRCON) {
				*hnext = SAFETY*h*pow(errmax, PGROW);
			} else {
				*hnext = 5.0*h;
			}
			*x += h;
			y[0] = ytemp[0];
			y[1] = ytemp[1];
			y[2] = ytemp[2];
			break; //return
		}
	}
}

inline void steeper_rkck(double *xpoint, double *grdt, double h0, double *xout, double *xerr){

  static const double b21 = 1.0/5.0;
  static const double b31 = 3.0/40.0;
  static const double b32 = 9.0/40.0;
  static const double b41 = 3.0/10.0;
  static const double b42 = -9.0/10.0;
  static const double b43 = 6.0/5.0;
  static const double b51 = -11.0/54.0;
  static const double b52 = 5.0/2.0;
  static const double b53 = -70.0/27.0;
  static const double b54 = 35.0/27.0;
  static const double b61 = 1631.0/55296.0;
  static const double b62 = 175.0/512.0;
  static const double b63 = 575.0/13824.0;
  static const double b64 = 44275.0/110592.0;
  static const double b65 = 253.0/4096.0;
  static const double c1 = 37.0/378.0;
  static const double c3 = 250.0/621.0;
  static const double c4 = 125.0/594.0;
  static const double c6 = 512.0/1771.0;
  double dc1 = c1-(2825.0/27648.0);
  double dc3 = c3-(18575.0/48384.0);
  double dc4 = c4-(13525.0/55296.0);
  double dc5 = -277.0/14336.0;
  double dc6 = c6-(1.0/4.0);
  
  double rho, grad[3], gradmod;

  xout[0] = xpoint[0] + h0*b21*grdt[0];
  xout[1] = xpoint[1] + h0*b21*grdt[1];
  xout[2] = xpoint[2] + h0*b21*grdt[2];

  double ak2[3] = {0.0};
  rho_grad(xout, &rho, grad, &gradmod);
  ak2[0] = grad[0];
  ak2[1] = grad[1];
  ak2[2] = grad[2];
  xout[0] = xpoint[0] + h0*(b31*grdt[0]+b32*ak2[0]);
  xout[1] = xpoint[1] + h0*(b31*grdt[1]+b32*ak2[1]);
  xout[2] = xpoint[2] + h0*(b31*grdt[2]+b32*ak2[2]);
  
  double ak3[3] = {0.0};
  rho_grad(xout, &rho, grad, &gradmod);
  ak3[0] = grad[0]; 
  ak3[1] = grad[1]; 
  ak3[2] = grad[2]; 
  xout[0] = xpoint[0] + h0*(b41*grdt[0]+b42*ak2[0]+b43*ak3[0]);
  xout[1] = xpoint[1] + h0*(b41*grdt[1]+b42*ak2[1]+b43*ak3[1]);
  xout[2] = xpoint[2] + h0*(b41*grdt[2]+b42*ak2[2]+b43*ak3[2]);
  
  double ak4[3] = {0.0};
  rho_grad(xout, &rho, grad, &gradmod);
  ak4[0] = grad[0];
  ak4[1] = grad[1];
  ak4[2] = grad[2];
  xout[0] = xpoint[0] + h0*(b51*grdt[0]+b52*ak2[0]+b53*ak3[0]+b54*ak4[0]);
  xout[1] = xpoint[1] + h0*(b51*grdt[1]+b52*ak2[1]+b53*ak3[1]+b54*ak4[1]);
  xout[2] = xpoint[2] + h0*(b51*grdt[2]+b52*ak2[2]+b53*ak3[2]+b54*ak4[2]);
  
  double ak5[3] = {0.0};
  rho_grad(xout, &rho, grad, &gradmod);
  ak5[0] = grad[0];
  ak5[1] = grad[1];
  ak5[2] = grad[2];
  xout[0] = xpoint[0] + h0*(b61*grdt[0]+b62*ak2[0]+b63*ak3[0]+b64*ak4[0]+b65*ak5[0]);
  xout[1] = xpoint[1] + h0*(b61*grdt[1]+b62*ak2[1]+b63*ak3[1]+b64*ak4[1]+b65*ak5[1]);
  xout[2] = xpoint[2] + h0*(b61*grdt[2]+b62*ak2[2]+b63*ak3[2]+b64*ak4[2]+b65*ak5[2]);
  
  double ak6[3] = {0.0};
  rho_grad(xout, &rho, grad, &gradmod);
  ak6[0] = grad[0];
  ak6[1] = grad[1];
  ak6[2] = grad[2];
  xout[0] = xpoint[0] + h0*(c1*grdt[0]+c3*ak3[0]+c4*ak4[0]+c6*ak6[0]);
  xout[1] = xpoint[1] + h0*(c1*grdt[1]+c3*ak3[1]+c4*ak4[1]+c6*ak6[1]);
  xout[2] = xpoint[2] + h0*(c1*grdt[2]+c3*ak3[2]+c4*ak4[2]+c6*ak6[2]);

  xerr[0] = h0*(dc1*grdt[0]+dc3*ak3[0]+dc4*ak4[0]+dc5*ak5[0]+dc6*ak6[0]);
  xerr[1] = h0*(dc1*grdt[1]+dc3*ak3[1]+dc4*ak4[1]+dc5*ak5[1]+dc6*ak6[1]);
  xerr[2] = h0*(dc1*grdt[2]+dc3*ak3[2]+dc4*ak4[2]+dc5*ak5[2]+dc6*ak6[2]);

	}


bool checkcp(const double *x, int *nuc){

  int i;
  bool iscp = false;
  double rho, grad[3], gradmod;

  *nuc = -2;
  rho_grad(x, &rho, grad, &gradmod);

  for (i=0; i<ncent_; i++){
    if (fabs(x[0]-xyzrho_[i*3+0]) < epsiscp_ &&
        fabs(x[1]-xyzrho_[i*3+1]) < epsiscp_ &&
        fabs(x[2]-xyzrho_[i*3+2]) < epsiscp_){
      iscp = true;
      *nuc = i;
      return iscp;
    }
  }

  // Put in the begining
  if (gradmod <= GRADEPS){
    iscp = true;
    if (rho <= RHOEPS){
      *nuc = -1;
    }
  } 

  return iscp;
  
}
 
inline void rho_grad(const double *p, double *rho, double *grad, double *gradmod){

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
  for (ic=0; ic<ncent_; ic++){
    xcoor[0] = p[0] - coords_[ic*3+0];
    xcoor[1] = p[1] - coords_[ic*3+1];
    xcoor[2] = p[2] - coords_[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
  	for (m=0; m<ngroup_[ic]; m++){
    	if (dis2 >= rcutte_[m*ncent_+ic]){continue;}
			int idx1 = m*(ngtoh*ncent_) + 0*ncent_ + ic;
      int k = nuexp_[idx1];
      double ori = -oexp_[k-1];
      double dp2 = ori + ori;
      double aexp = exp(ori*dis2);
      for (jj=0; jj<nzexp_[m*ncent_+ic]; jj++){
				int idx2 = m*(ngtoh*ncent_) + jj*ncent_ + ic;
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

void print_mole(){
	int i;
  for (i=0; i<ncent_; i++) {
	  printf("Coordinates for atom %d %f %f %f\n", i, 
    coords_[i*3+0], coords_[i*3+1], coords_[i*3+2]);
  }
}

// For a 3D matrix L by N by M:
// matrix[ i ][ j ][ k ] = array[ i*(N*M) + j*M + k ]
void print_basis(){
  int i,j,k;
  int suma = 0;
  printf("\n*Printing basis set info\n");
  printf("Number of orbitals %d\n", nmo_);
  printf("Number of primitives %d\n", nprims_);
  for (i=0; i<ncent_; i++) {
    suma = 0;
    printf("ngroup %d %d\n",i,ngroup_[i]);
    for (j=0; j<ngroup_[i]; j++) {
      printf("nzexp %d %d\n",j,nzexp_[j*ncent_+i]);
      printf("rcutte %d %f\n",j,rcutte_[j*ncent_+i]);
      suma += nzexp_[j*ncent_+i];
      for (k=0; k<nzexp_[j*ncent_+i]; k++) {
				int idx = j*(ngtoh*ncent_) + k*ncent_ + i;
        printf("nuexp, oexp %d %f\n", nuexp_[idx], oexp_[nuexp_[idx]-1]);
      }
    }
    printf("suma %d\n", suma);
  }
}

void init_nlm(){

  int i,j;

  for (i=0; i<56; i++) {
    for (j=0; j<3; j++) {
      nlm[i][j] = 0;
    }
  }

  // p's
  nlm[1][0] =  1;
  nlm[2][1] =  1;
  nlm[3][2] =  1;

  // d's
  nlm[4][0] =  2;
  nlm[5][1] =  2;
  nlm[6][2] =  2;
  nlm[7][0] =  1;
  nlm[7][1] =  1;
  nlm[8][0] =  1;
  nlm[8][2] =  1;
  nlm[9][1] =  1;
  nlm[9][2] =  1;

  // f's
  nlm[10][0] = 3;
  nlm[11][1] = 3;
  nlm[12][2] = 3;
  nlm[13][0] = 2;
  nlm[13][1] = 1;
  nlm[14][0] = 2;
  nlm[14][2] = 1;
  nlm[15][1] = 2;
  nlm[15][2] = 1;
  nlm[16][0] = 1;
  nlm[16][1] = 2;
  nlm[17][0] = 1;
  nlm[17][2] = 2;
  nlm[18][1] = 1;
  nlm[18][2] = 2;
  nlm[19][0] = 1;
  nlm[19][1] = 1;
  nlm[19][2] = 1;

  // g's
  nlm[20][0] = 4;
  nlm[21][1] = 4;
  nlm[22][2] = 4;
  nlm[23][0] = 3;
  nlm[23][1] = 1;
  nlm[24][0] = 3;
  nlm[24][2] = 1;
  nlm[25][0] = 1;
  nlm[25][1] = 3;
  nlm[26][1] = 3;
  nlm[26][2] = 1;
  nlm[27][0] = 1;
  nlm[27][2] = 3;
  nlm[28][1] = 1;
  nlm[28][2] = 3;
  nlm[29][0] = 2;
  nlm[29][1] = 2;
  nlm[30][0] = 2;
  nlm[30][2] = 2;
  nlm[31][1] = 2;
  nlm[31][2] = 2;
  nlm[32][0] = 2;
  nlm[32][1] = 1;
  nlm[32][2] = 1;
  nlm[33][0] = 1;
  nlm[33][1] = 2;
  nlm[33][2] = 1;
  nlm[34][0] = 1;
  nlm[34][1] = 1;
  nlm[34][2] = 2;

  // h's
  nlm[35][0] = 0;
  nlm[35][1] = 0;
  nlm[35][2] = 5;
  nlm[36][0] = 0;
  nlm[36][1] = 1;
  nlm[36][2] = 4;
  nlm[37][0] = 0;
  nlm[37][1] = 2;
  nlm[37][2] = 3;
  nlm[38][0] = 0;
  nlm[38][1] = 3;
  nlm[38][2] = 2;
  nlm[39][0] = 0;
  nlm[39][1] = 4;
  nlm[39][2] = 1;
  nlm[40][0] = 0;
  nlm[40][1] = 5;
  nlm[40][2] = 0;
  nlm[41][0] = 1;
  nlm[41][1] = 0;
  nlm[41][2] = 4;
  nlm[42][0] = 1;
  nlm[42][1] = 1;
  nlm[42][2] = 3;
  nlm[43][0] = 1;
  nlm[43][1] = 2;
  nlm[43][2] = 2;
  nlm[44][0] = 1;
  nlm[44][1] = 3;
  nlm[44][2] = 1;
  nlm[45][0] = 1;
  nlm[45][1] = 4;
  nlm[45][2] = 0;
  nlm[46][0] = 2;
  nlm[46][1] = 0;
  nlm[46][2] = 3;
  nlm[47][0] = 2;
  nlm[47][1] = 1;
  nlm[47][2] = 2;
  nlm[48][0] = 2;
  nlm[48][1] = 2;
  nlm[48][2] = 1;
  nlm[49][0] = 2;
  nlm[49][1] = 3;
  nlm[49][2] = 0;
  nlm[50][0] = 3;
  nlm[50][1] = 0;
  nlm[50][2] = 2;
  nlm[51][0] = 3;
  nlm[51][1] = 1;
  nlm[51][2] = 1;
  nlm[52][0] = 3;
  nlm[52][1] = 2;
  nlm[52][2] = 0;
  nlm[53][0] = 4;
  nlm[53][1] = 0;
  nlm[53][2] = 1;
  nlm[54][0] = 4;
  nlm[54][1] = 1;
  nlm[54][2] = 0;
  nlm[55][0] = 5;
  nlm[55][1] = 0;
  nlm[55][2] = 0;
}

inline void cerror(const char *text){
	fprintf(stderr,"Error %s\n", text);
	exit(1);
}

