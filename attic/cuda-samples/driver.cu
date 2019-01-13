
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "driver.h"

// Use shared memory and constant/restricted

__device__ __forceinline__ double rho(double *point, double *coeff,
                      double *occ, double *oxp, int *cen, int *typ,
                      double *xyz){

  double gun[20], trho;  // TODO: change gun
  double xcoor[3], fun[3];
  int i, j, it[3];

  trho = 0.0;
  for (i=0; i<dnmo_; i++) {
    gun[i] = 0.0;
  }

  // Simpler version over primitive, no screaning, no cuttof !!
  for (i=0; i<dnprim_; i++) {
    int ic = cen[i]-1;
    int itip = typ[i]-1;
    it[0] = dnlm_[itip][0];
    it[1] = dnlm_[itip][1];
    it[2] = dnlm_[itip][2];
    double ori = -oxp[i];
    xcoor[0] = point[0] - xyz[ic*3+0];
    xcoor[1] = point[1] - xyz[ic*3+1];
    xcoor[2] = point[2] - xyz[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
    double aexp = exp(ori*dis2);
    for (j=0; j<3; j++) {
      int n = it[j];
      double x = xcoor[j];
      if (n==0) {
        fun[j] = 1.0;
      }
      else if (n==1) {
        fun[j] = x;
      }
      else if (n==2) {
        fun[j] = x*x;
      }
      else if (n==3) {
        fun[j] = x*x*x;
      }
      else if (n==4) {
        fun[j] = x*x*x*x;
      }
      else if (n==5) {
        fun[j] = x*x*x*x*x;
      }
    }
    double f123 = fun[0]*fun[1]*fun[2]*aexp;
    for (j=0; j<dnmo_; j++) {
      double cfj = coeff[i*dnmo_+j];
      gun[j] += cfj*f123;
    }
  }

  // Run again over orbitals
  for (i=0; i<dnmo_; i++) {
    trho += occ[i]*gun[i]*gun[i];
  }

  return trho;

}


__global__ void compute(double *out, double *points, double *coeff,
                        double *occ, double *oxp, int *cen, int *typ,
                        double *xyz, const int npoints){

  const unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
  //const unsigned int numThreads = blockDim.x * gridDim.x;

	/*int i, j;
	__shared__ double coeffs[1200];
  __shared__ double occs[200];
  __shared__ double oxps[200];
	__shared__ int cens[200];
	__shared__ int typs[200];
	__shared__ double xyzs[200];

  for (i=0; i<dnprim_; i++) {
		cens[i] = cen[i];
		typs[i] = typ[i];
		oxps[i] = oxp[i];
	}
  for (i=0; i<dnmo_; i++) {
		occs[i] = occ[i];
  	for (j=0; j<dnprim_; j++) {
    	coeffs[j*dnmo_+i] = coeff[j*dnmo_+i];
		}
	}
  for (i=0; i<dnatm_; i++) {
    xyzs[i*3+0] = xyz[i*3+0];
    xyzs[i*3+1] = xyz[i*3+1];
    xyzs[i*3+2] = xyz[i*3+2];
  }*/

  if (idx < npoints) {
    double point[3];
    point[0] = points[idx*3+0];
    point[1] = points[idx*3+1];
    point[2] = points[idx*3+2];
    double trho = rho(point, coeff, occ, oxp, cen, typ, xyz);
    out[idx] = trho;
  }

  //__syncthreads();

}

//for a 3D matrix L by N by M:
//matrix[ i ][ j ][ k ] = array[ i*(N*M) + j*M + k ]
extern "C" {
  void basis(const int nmo, const int nprim, const int natm,
             const int *ngroup, const int *nzexp, const int *nuexp,
             const double *rcutte){

    const int mgrp = 200;
    const int ngtoh = 21;
    int i,j,k;
    int suma = 0;

    for (i=0; i<natm; i++) {
      suma = 0;
      printf("ngroup %d %d\n",i,ngroup[i]);
      for (j=0; j<ngroup[i]; j++) {
        printf("nzexp %d %d\n",j,nzexp[j*natm+i]);
        printf("Rcutte %d %f\n",j,rcutte[j*natm+i]);
        suma += nzexp[j*natm+i];
        for (k=0; k<nzexp[j*natm+i]; k++) {
					int idx = j*(ngtoh*natm) + k*natm + i;
          printf("nuexp %d\n", nuexp[idx]);
        }
      }
      printf("suma %d\n", suma);
    }
	
	}
}


extern "C" {
  void driver(const int nmo, const int nprim, const int natm,
              const int npoints, const int *icen, const int *ityp,
              const double *oexp, const double *mo_coeff,
              const double *mo_occ, const double *coords,
              const double *bcoords, double *brho) {

    int i, j;
    const int blockSize = 1024;
    int gridSize;

    init_nlm();

    // Basis info
    natm_ = natm;
    nprim_ = nprim;
    nmo_ = nmo;
    mo_coeff_ = (double *) malloc(sizeof(double)*nmo_*nprim_);
    assert(mo_coeff_ != NULL);
    mo_occ_ = (double *) malloc(sizeof(double)*nmo_);
    assert(mo_occ_ != NULL);
    icen_ = (int *) malloc(sizeof(int)*nprim_);
    assert(icen_ != NULL);
    ityp_ = (int *) malloc(sizeof(int)*nprim_);
    assert(ityp_ != NULL);
    oexp_ = (double *) malloc(sizeof(double)*nprim_);
    assert(oexp_ != NULL);
    for (i=0; i<nprim_; i++) {
      icen_[i] = icen[i];
      ityp_[i] = ityp[i];
      oexp_[i] = oexp[i];
    }
    for (i=0; i<nmo_; i++) {                      // Orbital
      mo_occ_[i] = mo_occ[i];
      for (j=0; j<nprim_; j++) {
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

    //double point[3];
    //point[0] = 0.0;
    //point[1] = 0.0;
    //point[2] = 0.0;
    //double rho = point_rho(point);
    //printf("The value of rho is %f\n", rho);

    // GPU basis info
    cudaMalloc((void **)&dbrho, npoints*sizeof(double));
    cudaMalloc((void **)&dbcoords, 3*npoints*sizeof(double));

    cudaMalloc((void **)&dmo_coeff_, nmo_*nprim_*sizeof(double));
    cudaMalloc((void **)&dmo_occ_, nmo_*sizeof(double));
    cudaMalloc((void **)&dicen_, nprim_*sizeof(int));
    cudaMalloc((void **)&dityp_, nprim_*sizeof(int));
    cudaMalloc((void **)&doexp_, nprim_*sizeof(double));
    cudaMalloc((void **)&dcoords_, 3*natm_*sizeof(double));

    cudaMemcpy(dbcoords, bcoords, npoints*3*sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(dmo_coeff_, mo_coeff_, nmo_*nprim_*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dmo_occ_, mo_occ_, nmo_*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dicen_, icen_, nprim_*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dityp_, ityp_, nprim_*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(doexp_, oexp_, nprim_*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dcoords_, coords_, 3*natm_*sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(dnprim_, &nprim_, sizeof(int));
    cudaMemcpyToSymbol(dnmo_, &nmo_, sizeof(int));
    cudaMemcpyToSymbol(dnatm_, &natm_, sizeof(int));
    cudaMemcpyToSymbol(dnlm_, &nlm, 56*3*sizeof(int));

    gridSize = (npoints + blockSize - 1) / blockSize;

    float timer;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    compute<<<gridSize,blockSize>>>(dbrho, dbcoords, dmo_coeff_, dmo_occ_, doexp_, dicen_, dityp_, dcoords_, npoints);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&timer, start, stop);
    printf("GPU elapsed time:  %3.3f s \n", timer/1000);

    cudaMemcpy(brho, dbrho, npoints*sizeof(double), cudaMemcpyDeviceToHost);
    //for (i = 0; i<2; ++i){
    //  printf("%f\n", brho[i]);
    //}

    free(mo_coeff_);
    free(mo_occ_);
    free(oexp_);
    free(coords_);
    free(icen_);
    free(ityp_);
    cudaFree(dmo_coeff_);
    cudaFree(dmo_occ_);
    cudaFree(doexp_);
    cudaFree(dcoords_);
    cudaFree(dicen_);
    cudaFree(dityp_);
    cudaFree(dbrho);
    cudaFree(dbcoords);

  }
}


// Slow versions no shells
inline double point_rho(double *p)
{

  double gun[nmo_], rho;
  double xcoor[3], fun[3];
  int i, j, it[3];

  rho = 0.0;
  for (i=0; i<nmo_; i++) {
    gun[i] = 0.0;
  }

  // Simpler version over primitive, no screaning, no cuttof !!
  for (i=0; i<nprim_; i++) {
    int ic = icen_[i]-1;
    int itip = ityp_[i]-1;
    it[0] = nlm[itip][0];
    it[1] = nlm[itip][1];
    it[2] = nlm[itip][2];
    double ori = -oexp_[i];
    xcoor[0] = p[0] - coords_[ic*3+0];
    xcoor[1] = p[1] - coords_[ic*3+1];
    xcoor[2] = p[2] - coords_[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
    double aexp = exp(ori*dis2);
    for (j=0; j<3; j++) {
      int n = it[j];
      double x = xcoor[j];
      if (n==0) {
        fun[j] = 1.0;
      }
      else if (n==1) {
        fun[j] = x;
      }
      else if (n==2) {
        fun[j] = x*x;
      }
      else if (n==3) {
        fun[j] = x*x*x;
      }
      else if (n==4) {
        fun[j] = x*x*x*x;
      }
      else if (n==5) {
        fun[j] = x*x*x*x*x;
      }
    }
    double f123 = fun[0]*fun[1]*fun[2]*aexp;
    for (j=0; j<nmo_; j++) {
      double cfj = mo_coeff_[i*nmo_+j];
      gun[j] += cfj*f123;
    }
  }

  // Run again over orbitals
  for (i=0; i<nmo_; i++) {
    rho += mo_occ_[i]*gun[i]*gun[i];
  }

  return rho;

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

