#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "gsurf.h"

extern "C" { 
void gsurf_driver(const int nmo, 
                  const int nprims,
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

  const int blockSize = 128;
  const int gridSize = (npang + blockSize - 1) / blockSize; 

  init_nlm();

  cudaMemcpyToSymbol(dnprims_, &nprims, sizeof(int));
  cudaMemcpyToSymbol(dnmo_, &nmo, sizeof(int));
  cudaMemcpyToSymbol(dncent_, &natm, sizeof(int));
  cudaMemcpyToSymbol(dnlm, &nlm, 56*3*sizeof(int));
  cudaMemcpyToSymbol(dnpang_, &npang, sizeof(int));
  cudaMemcpyToSymbol(dinuc_, &inuc, sizeof(int));
  cudaMemcpyToSymbol(depsiscp_, &epsiscp, sizeof(double));
  cudaMemcpyToSymbol(dntrial_, &ntrial, sizeof(int));
  cudaMemcpyToSymbol(depsroot_, &epsroot, sizeof(double));
  cudaMemcpyToSymbol(depsilon_, &epsilon, sizeof(double));
  cudaMemcpyToSymbol(drmaxsurf_, &rmaxsurf, sizeof(double));
  cudaMemcpyToSymbol(dbackend_, &backend, sizeof(int));
  cudaMemcpyToSymbol(dstep_, &step, sizeof(double));
  cudaMemcpyToSymbol(dmstep_, &mstep, sizeof(int));
  cudaMemcpyToSymbol(dsteeper_, &backend, sizeof(int));
	double xnuc[3];
  xnuc[0] = xyzrho[inuc*3+0];
  xnuc[1] = xyzrho[inuc*3+1]; 
  xnuc[2] = xyzrho[inuc*3+2];
  cudaMemcpyToSymbol(dxnuc_, &xnuc, 3*sizeof(double));

  cudaMalloc((void **)&dmo_coeff_, nmo*nprims*sizeof(double));
  cudaMalloc((void **)&dmo_occ_, nmo*sizeof(double));
  cudaMalloc((void **)&dityp_, nprims*sizeof(int));
  cudaMalloc((void **)&doexp_, nprims*sizeof(double));
  cudaMalloc((void **)&dcoords_, 3*natm*sizeof(double));
  cudaMalloc((void **)&dxyzrho_, 3*natm*sizeof(double));
  cudaMalloc((void **)&dngroup_, natm*sizeof(int));
  cudaMalloc((void **)&dnzexp_, natm*mgrp*sizeof(int));
  cudaMalloc((void **)&dnuexp_, natm*mgrp*ngtoh*sizeof(int));
  cudaMalloc((void **)&drcutte_, natm*mgrp*sizeof(double));
  cudaMalloc((void **)&drpru_, ntrial*sizeof(double));
  cudaMalloc((void **)&dct_, npang*sizeof(double));
  cudaMalloc((void **)&dst_, npang*sizeof(double));
  cudaMalloc((void **)&dcp_, npang*sizeof(double));
  cudaMalloc((void **)&dsp_, npang*sizeof(double));
  cudaMalloc((void **)&dnlimsurf_, npang*sizeof(int));
  cudaMalloc((void **)&drsurf_, ntrial*npang*sizeof(double));

  cudaMemcpy(dmo_coeff_, mo_coeff, nmo*nprims*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dmo_occ_, mo_occ, nmo*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dityp_, ityp, nprims*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(doexp_, oexp, nprims*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dcoords_, coords, 3*natm*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dxyzrho_, xyzrho, 3*natm*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dngroup_, ngroup, natm*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dnzexp_, nzexp, natm*mgrp*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dnuexp_, nuexp, natm*mgrp*ngtoh*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(drcutte_, rcutte, natm*mgrp*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(drpru_, rpru, ntrial*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dct_, ct, npang*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dst_, st, npang*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dcp_, cp, npang*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dsp_, sp, npang*sizeof(double), cudaMemcpyHostToDevice);

	float timer;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  printf("GPU dim grid,block %d %d\n", gridSize, blockSize);
  compute<<<gridSize,blockSize>>>(dityp_,doexp_,dngroup_,dnzexp_,dnuexp_,
                                  drcutte_,dmo_coeff_,dmo_occ_,dcoords_,
                                  dxyzrho_,drpru_,dct_,dst_,dcp_,dsp_,
                                  dnlimsurf_,drsurf_);
  cudaError cuerr = cudaGetLastError();
 	printf("CUDA Error : %s\n", cudaGetErrorString(cuerr)); 
  cudaDeviceSynchronize();
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timer, start, stop);
  cudaMemcpy(nlimsurf, dnlimsurf_, npang*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(rlimsurf, drsurf_, ntrial*npang*sizeof(double), cudaMemcpyDeviceToHost);
	//for (i=0; i<npang_; i++){
  //  nlimsurf[i] = nlimsurf_[i];
	//  for (j=0; j<ntrial_; j++){
  //    rlimsurf[i*ntrial_+j] = rsurf_[i*ntrial_+j];
  //  }
  //}
  printf("GPU elapsed time: %3.3f s\n", timer/1000);

  cudaFree(dngroup_);
  cudaFree(dnuexp_);
  cudaFree(dnzexp_);
  cudaFree(drcutte_);
  cudaFree(dmo_coeff_);
  cudaFree(dmo_occ_);
  cudaFree(doexp_);
  cudaFree(dcoords_);
  cudaFree(dityp_);
  cudaFree(dcoords_);
  cudaFree(dxyzrho_);
  cudaFree(drcutte_);
  cudaFree(dcp_);
  cudaFree(dsp_);
  cudaFree(dct_);
  cudaFree(dst_);


}}

__device__ bool __forceinline__ checkcp(double *x, int *nuc,
                        const int *ityp,            
                        const double *oexp,         
                        const int *ngroup,          
                        const int *nzexp,           
                        const int *nuexp,           
                        const double *rcutte,       
                        const double *mo_coeff,     
                        const double *mo_occ,       
                        const double *coords,
            						const double *xyzrho){

  int i;
  bool iscp = false;
  double rho, grad[3], gradmod;

  *nuc = -2;
  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, x, &rho, grad, &gradmod);

  for (i=0; i<dncent_; i++){
    if (fabs(x[0]-xyzrho[i*3+0]) < depsiscp_ &&
        fabs(x[1]-xyzrho[i*3+1]) < depsiscp_ &&
        fabs(x[2]-xyzrho[i*3+2]) < depsiscp_){
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

__global__ void compute(const int *ityp,         
                        const double *oexp,      
                        const int *ngroup,       
                        const int *nzexp,        
                        const int *nuexp,        
                        const double *rcutte,    
                        const double *mo_coeff,  
                        const double *mo_occ,    
                        const double *coords,    
                        const double *xyzrho,    
                        const double *rpru,
				                const double *ct,         
				                const double *st,         
				                const double *cp,         
				                const double *sp,          
				                int *nlimsurf,         
				                double *rsurf){         
                        
  const unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < dnpang_) {
  	int j;
  	double xsurf[11][3]; //dntrial=11 
	  int isurf[11][2];
	  int nintersec = 0;
	  double xpoint[3];
  	if (dncent_ == 1){  
      nlimsurf[idx] = 1;
      int offset = idx*dntrial_;
      rsurf[offset] = drmaxsurf_;
      return;
    }
    nintersec = 0;
    double cost = ct[idx];
    double sintcosp = st[idx]*cp[idx];
    double sintsinp = st[idx]*sp[idx];
    int ia = dinuc_, ib;
    double ra = 0.0, rb;
	  for (j=0; j<dntrial_; j++){
      double ract = rpru[j];
      xpoint[0] = dxnuc_[0] + ract*sintcosp; 
      xpoint[1] = dxnuc_[1] + ract*sintsinp; 
      xpoint[2] = dxnuc_[2] + ract*cost;     
      //TODO: Better Check for error
      int ier = odeint(xpoint, dstep_, depsilon_, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords,xyzrho); 
      if (ier == 1) {
        printf("too short steep on odeint\n");
      } else if (ier == 2) {
        printf("too many steeps on odeint\n");
			} else if (ier == 4) {
        printf("nna on odeint\n");
			}
      bool good = checkcp(xpoint, &ib, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords,xyzrho);  
      rb = ract;
      if (ib != ia && (ia == dinuc_ || ib == dinuc_)){
        if (ia != dinuc_ || ib != -1){
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
      xin[0] = dxnuc_[0] + ra*sintcosp;
      xin[1] = dxnuc_[1] + ra*sintsinp;
      xin[2] = dxnuc_[2] + ra*cost;
      xfin[0] = dxnuc_[0] + rb*sintcosp;
      xfin[1] = dxnuc_[1] + rb*sintsinp;
      xfin[2] = dxnuc_[2] + rb*cost;
      while (fabs(ra-rb) > depsroot_){
        double xmed[3];
        xmed[0] = 0.5*(xfin[0] + xin[0]);   
        xmed[1] = 0.5*(xfin[1] + xin[1]);   
        xmed[2] = 0.5*(xfin[2] + xin[2]);   
        double rm = 0.5*(ra + rb);
        xpoint[0] = xmed[0];
        xpoint[1] = xmed[1];
        xpoint[2] = xmed[2];
        int im;
        int ier = odeint(xpoint, dstep_, depsilon_, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords,xyzrho);  
      	if (ier == 1) {
        	printf("too short steep on odeint\n");
	      } else if (ier == 2) {
        	printf("too many steeps on odeint\n");
				} else if (ier == 4) {
        	printf("nna on odeint\n");
				}
        bool good = checkcp(xpoint, &im, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords,xyzrho);   
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
          if (ia == dinuc_){
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
    nlimsurf[idx] = nintersec; 
	  for (j=0; j<nintersec; j++){
      int offset = idx*dntrial_+j;
      rsurf[offset] = xsurf[j][2];
    }
    if (nintersec%2 == 0){
      nintersec += 1;
      nlimsurf[idx] = nintersec;
      int offset = idx*dntrial_+(nintersec-1);
      rsurf[offset] = drmaxsurf_;
    }
    //printf("#* %d %d %.6f %.6f %.6f %.6f ",i,nlimsurf_[i],ct_[i],st_[i],cp_[i],sp_[i]);
	  //for (j=0; j<nlimsurf_[i]; j++){
    // printf(" %.6f ",rsurf_[i*ntrial_+j]);
    //}
    //printf("\n");
	}
}

//ier = 0 (correct), 1 (short step), 2 (too many iterations), 
//      3 (infty), 4 (nna), 5(undef)
__device__ int odeint(double *ystart, double h1, double eps,
                      const int *ityp,            
                      const double *oexp,         
                      const int *ngroup,          
                      const int *nzexp,           
                      const int *nuexp,           
                      const double *rcutte,       
                      const double *mo_coeff,     
                      const double *mo_occ,       
                      const double *coords,
            					const double *xyzrho){

	double rho, gradmod, hnext;
  double grad[3] = {0.0};
	int ier = 0, i, nuc;
  double dydx[3], y[3], yscal[3];

	rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, ystart, &rho, grad, &gradmod);
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

	for (i=0; i<dmstep_; i++){
	  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, y, &rho, dydx, &gradmod);
		yscal[0] = fmax(fabs(y[0]) + fabs(dydx[0]*h) + TINY, eps);
		yscal[1] = fmax(fabs(y[1]) + fabs(dydx[1]*h) + TINY, eps);
		yscal[2] = fmax(fabs(y[2]) + fabs(dydx[2]*h) + TINY, eps);
		if ((x+h-x2)*(x+h-x1) > 0.0) h = x2 - x;
		rkqs(y, dydx, &x, h, eps, yscal, &hnext, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords);
		if ((x-x2)*(x2-x1) >= 0.0 || checkcp(y, &nuc, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords,xyzrho)){
			ystart[0] = y[0];
			ystart[1] = y[1];
			ystart[2] = y[2];
      ier = 0;
			return ier;
		}
		if (fabs(hnext) <= hmin) printf("Step size too small in odeint\n");
		if (i == (dmstep_-1)) printf("Reached max steps in odeint\n");
		h = hnext;
  }
    
	// Test if the point is far from RMAXSURF from current atom. 
  double a1 = y[0] - dxnuc_[0];
  double a2 = y[1] - dxnuc_[1];
  double a3 = y[2] - dxnuc_[2];
  if ((a1*a1+a2*a2+a3*a3)>=5.0*5.0){
    ier = 3;
  } else { 
  	printf("NNA at %f %f %f\n", y[0], y[1], y[2]);
  }

  return ier;
}

__device__ void rkqs(double *y, double *dydx, double *x, 
                     double htry, double eps,
	                   double *yscal, double *hnext,
                     const int *ityp,            
                     const double *oexp,         
                     const int *ngroup,          
                     const int *nzexp,           
                     const int *nuexp,           
                     const double *rcutte,       
                     const double *mo_coeff,     
                     const double *mo_occ,       
                     const double *coords){
	
  int i;
	double yerr[3], ytemp[3], errmax, xnew, htemp;

	double h = htry;
  *hnext = 0.0;
  errmax = 0.0;

	for (;;){
		steeper_rkck(y, dydx, h, ytemp, yerr, ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords);
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
        printf("stepsize underflow in rkqs\n");
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

__device__ __forceinline__ void steeper_rkck(double *xpoint, double *grdt, double h0, double *xout, double *xerr,
                                             const int *ityp,            
                                             const double *oexp,         
                                             const int *ngroup,          
                                             const int *nzexp,           
                                             const int *nuexp,           
                                             const double *rcutte,       
                                             const double *mo_coeff,     
                                             const double *mo_occ,       
                                             const double *coords){

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
  static double dc1 = c1-(2825.0/27648.0);
  static double dc3 = c3-(18575.0/48384.0);
  static double dc4 = c4-(13525.0/55296.0);
  static double dc5 = -277.0/14336.0;
  static double dc6 = c6-(1.0/4.0);
  
  double rho, grad[3], gradmod;

  xout[0] = xpoint[0] + h0*b21*grdt[0];
  xout[1] = xpoint[1] + h0*b21*grdt[1];
  xout[2] = xpoint[2] + h0*b21*grdt[2];

  double ak2[3] = {0.0};
  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, xout, &rho, grad, &gradmod);
  ak2[0] = grad[0];
  ak2[1] = grad[1];
  ak2[2] = grad[2];
  xout[0] = xpoint[0] + h0*(b31*grdt[0]+b32*ak2[0]);
  xout[1] = xpoint[1] + h0*(b31*grdt[1]+b32*ak2[1]);
  xout[2] = xpoint[2] + h0*(b31*grdt[2]+b32*ak2[2]);
  
  double ak3[3] = {0.0};
  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, xout, &rho, grad, &gradmod);
  ak3[0] = grad[0]; 
  ak3[1] = grad[1]; 
  ak3[2] = grad[2]; 
  xout[0] = xpoint[0] + h0*(b41*grdt[0]+b42*ak2[0]+b43*ak3[0]);
  xout[1] = xpoint[1] + h0*(b41*grdt[1]+b42*ak2[1]+b43*ak3[1]);
  xout[2] = xpoint[2] + h0*(b41*grdt[2]+b42*ak2[2]+b43*ak3[2]);
  
  double ak4[3] = {0.0};
  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, xout, &rho, grad, &gradmod);
  ak4[0] = grad[0];
  ak4[1] = grad[1];
  ak4[2] = grad[2];
  xout[0] = xpoint[0] + h0*(b51*grdt[0]+b52*ak2[0]+b53*ak3[0]+b54*ak4[0]);
  xout[1] = xpoint[1] + h0*(b51*grdt[1]+b52*ak2[1]+b53*ak3[1]+b54*ak4[1]);
  xout[2] = xpoint[2] + h0*(b51*grdt[2]+b52*ak2[2]+b53*ak3[2]+b54*ak4[2]);
  
  double ak5[3] = {0.0};
  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, xout, &rho, grad, &gradmod);
  ak5[0] = grad[0];
  ak5[1] = grad[1];
  ak5[2] = grad[2];
  xout[0] = xpoint[0] + h0*(b61*grdt[0]+b62*ak2[0]+b63*ak3[0]+b64*ak4[0]+b65*ak5[0]);
  xout[1] = xpoint[1] + h0*(b61*grdt[1]+b62*ak2[1]+b63*ak3[1]+b64*ak4[1]+b65*ak5[1]);
  xout[2] = xpoint[2] + h0*(b61*grdt[2]+b62*ak2[2]+b63*ak3[2]+b64*ak4[2]+b65*ak5[2]);
  
  double ak6[3] = {0.0};
  rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,mo_occ,coords, xout, &rho, grad, &gradmod);
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

__device__ __forceinline__ void rho_grad(const int *ityp,        
                                         const double *oexp,     
                                         const int *ngroup,      
                                         const int *nzexp,       
                                         const int *nuexp,       
                                         const double *rcutte,   
                                         const double *mo_coeff, 
                                         const double *mo_occ,   
                                         const double *coords,   
                                         const double *point,    
                                         double *rho,            
                                         double *grad,           
                                         double *gradmod){       
                                         
  double gun[20];
	double gun1x[20], gun1y[20], gun1z[20]; // Change to nmo value
  double xcoor[3], fun[3], fun1[3];
  unsigned int ii, ic, m, jj, j, it[3];

  *rho = 0.0;
  *gradmod = 0.0;
  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = 0.0;
  for (ii=0; ii<dnmo_; ii++) {
    gun[ii] = 0.0;
    gun1x[ii] = 0.0;
    gun1y[ii] = 0.0;
    gun1z[ii] = 0.0;
  }
    
	// TODO: avoid natm loop and loop only over total shells
  for (ic=0; ic<dncent_; ic++){
    xcoor[0] = point[0] - coords[ic*3+0];
    xcoor[1] = point[1] - coords[ic*3+1];
    xcoor[2] = point[2] - coords[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
  	for (m=0; m<ngroup[ic]; m++){
    	if (dis2 >= rcutte[m*dncent_+ic]){continue;}
			unsigned int idx1 = m*(ngtoh*dncent_) + 0*dncent_ + ic;
      unsigned int k = nuexp[idx1];
      double ori = -oexp[k-1];
      double dp2 = ori + ori;
      double aexp = exp(ori*dis2);
      for (jj=0; jj<nzexp[m*dncent_+ic]; jj++){
				unsigned int idx2 = m*(ngtoh*dncent_) + jj*dncent_ + ic;
      	unsigned int i = nuexp[idx2];
        unsigned int itip = ityp[i-1]-1;
    		it[0] = dnlm[itip][0];
		    it[1] = dnlm[itip][1];
		    it[2] = dnlm[itip][2];
    		for (j=0; j<3; j++) {
		      unsigned int n = it[j];
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
    		for (j=0; j<dnmo_; j++) {
      		double cfj = mo_coeff[(i-1)*dnmo_+j];
          gun[j] += cfj*f123;
          gun1x[j] += cfj*fa;
          gun1y[j] += cfj*fb;
          gun1z[j] += cfj*fc;
				}
      }
    }
  }
  // Run again over orbitals
  for (ii=0; ii<dnmo_; ii++) {
    *rho += mo_occ[ii]*gun[ii]*gun[ii];
    grad[0] += mo_occ[ii]*gun[ii]*gun1x[ii];
    grad[1] += mo_occ[ii]*gun[ii]*gun1y[ii];
    grad[2] += mo_occ[ii]*gun[ii]*gun1z[ii];
  }
  grad[0] += grad[0];
  grad[1] += grad[1];
  grad[2] += grad[2];
  *gradmod = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
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

