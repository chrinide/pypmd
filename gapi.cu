#include <stdio.h> 
#include <math.h>

#include "gapi.h"

extern "C" {
void eval_gpu(const int nmo, 
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
              const int npoints,
              const double *points,
              double *output){

  const int blockSize = 1024;
  const int gridSize = (npoints + blockSize - 1) / blockSize; 

  init_nlm();

  cudaMemcpyToSymbol(dnprims, &nprims, sizeof(int));
  cudaMemcpyToSymbol(dnmo, &nmo, sizeof(int));
  cudaMemcpyToSymbol(dncent, &natm, sizeof(int));
  cudaMemcpyToSymbol(dnlm, &nlm, 56*3*sizeof(int));
  cudaMemcpyToSymbol(dnpoints, &npoints, sizeof(int));

  double *dmo_coeff, *dmo_occ, *doexp, *dcoords, *drcutte; 
  int *dityp, *dngroup, *dnzexp, *dnuexp;
	double *dpoints, *doutput;
  cudaMalloc((void **)&dmo_coeff, nmo*nprims*sizeof(double));
  cudaMalloc((void **)&dmo_occ, nmo*sizeof(double));
  cudaMalloc((void **)&dityp, nprims*sizeof(int));
  cudaMalloc((void **)&doexp, nprims*sizeof(double));
  cudaMalloc((void **)&dcoords, 3*natm*sizeof(double));
  cudaMalloc((void **)&dpoints, 3*npoints*sizeof(double));
  cudaMalloc((void **)&doutput, npoints*sizeof(double));
  cudaMalloc((void **)&dngroup, natm*sizeof(int));
  cudaMalloc((void **)&dnzexp, natm*mgrp*sizeof(int));
  cudaMalloc((void **)&dnuexp, natm*mgrp*ngtoh*sizeof(int));
  cudaMalloc((void **)&drcutte, natm*mgrp*sizeof(double));

  cudaMemcpy(dmo_coeff, mo_coeff, nmo*nprims*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dmo_occ, mo_occ, nmo*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dityp, ityp, nprims*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(doexp, oexp, nprims*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dcoords, coords, 3*natm*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dpoints, points, 3*npoints*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dngroup, ngroup, natm*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dnzexp, nzexp, natm*mgrp*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dnuexp, nuexp, natm*mgrp*ngtoh*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(drcutte, rcutte, natm*mgrp*sizeof(double), cudaMemcpyHostToDevice);
    
	float timer;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  printf("GPU dim grid,block %d %d\n", gridSize, blockSize);
  compute<<<gridSize,blockSize>>>(dityp,doexp,dngroup,dnzexp,dnuexp,drcutte,dmo_coeff,dmo_occ,dcoords,
                                  dpoints,doutput);
  cudaError cuerr = cudaGetLastError();
 	fprintf(stderr, "CUDA Error : %s\n", cudaGetErrorString(cuerr)); 
  cudaDeviceSynchronize();
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timer, start, stop);
  cudaMemcpy(output, doutput, npoints*sizeof(double), cudaMemcpyDeviceToHost);
  printf("GPU elapsed time: %3.3f s\n", timer/1000);

  cudaFree(dngroup);
  cudaFree(dnuexp);
  cudaFree(dnzexp);
  cudaFree(drcutte);
  cudaFree(dmo_coeff);
  cudaFree(dmo_occ);
  cudaFree(doexp);
  cudaFree(dcoords);
  cudaFree(dityp);
  cudaFree(dcoords);
  cudaFree(dpoints);
  cudaFree(doutput);

}}

__global__ void compute(const int *ityp,         
                        const double *oexp,      
                        const int *ngroup,       
                        const int *nzexp,        
                        const int *nuexp,        
                        const double *rcutte,    
                        const double *mo_coeff,  
                        const double *mo_occ,    
                        const double *coords,    
                        const double *points,    
                        double *output){         
                        
  const unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < dnpoints) {
    double point[3], grad[3], rho, gradmod;
    point[0] = points[idx*3+0];
    point[1] = points[idx*3+1];
    point[2] = points[idx*3+2];
    rho_grad(ityp,oexp,ngroup,nzexp,nuexp,rcutte,mo_coeff,
             mo_occ,coords,point,&rho,grad,&gradmod);
    output[idx] = rho;
  }
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
                                         
  double gun[6], gun1x[6], gun1y[6], gun1z[6]; // Change to nmo value
  double xcoor[3], fun[3], fun1[3];
  unsigned int ii, ic, m, jj, j, it[3];

  *rho = 0.0;
  //*gradmod = 0.0;
  //grad[0] = 0.0;
  //grad[1] = 0.0;
  //grad[2] = 0.0;
  for (ii=0; ii<dnmo; ii++) {
    gun[ii] = 0.0;
    //gun1x[ii] = 0.0;
    //gun1y[ii] = 0.0;
    //gun1z[ii] = 0.0;
  }
    
	// TODO: avoid natm loop and loop only over total shells
  for (ic=0; ic<dncent; ic++){
    xcoor[0] = point[0] - coords[ic*3+0];
    xcoor[1] = point[1] - coords[ic*3+1];
    xcoor[2] = point[2] - coords[ic*3+2];
    double dis2 = xcoor[0]*xcoor[0] + xcoor[1]*xcoor[1] + xcoor[2]*xcoor[2];
  	for (m=0; m<ngroup[ic]; m++){
    	if (dis2 >= rcutte[m*dncent+ic]){continue;}
			unsigned int idx1 = m*(ngtoh*dncent) + 0*dncent + ic;
      unsigned int k = nuexp[idx1];
      double ori = -oexp[k-1];
      double dp2 = ori + ori;
      double aexp = exp(ori*dis2);
      for (jj=0; jj<nzexp[m*dncent+ic]; jj++){
				unsigned int idx2 = m*(ngtoh*dncent) + jj*dncent + ic;
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
    		for (j=0; j<dnmo; j++) {
      		double cfj = mo_coeff[(i-1)*dnmo+j];
          gun[j] += cfj*f123;
          //gun1x[j] += cfj*fa;
          //gun1y[j] += cfj*fb;
          //gun1z[j] += cfj*fc;
				}
      }
    }
  }
  // Run again over orbitals
  for (ii=0; ii<dnmo; ii++) {
    *rho += mo_occ[ii]*gun[ii]*gun[ii];
    //grad[0] += mo_occ[ii]*gun[ii]*gun1x[ii];
    //grad[1] += mo_occ[ii]*gun[ii]*gun1y[ii];
    //grad[2] += mo_occ[ii]*gun[ii]*gun1z[ii];
  }
  //grad[0] += grad[0];
  //grad[1] += grad[1];
  //grad[2] += grad[2];
  //*gradmod = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
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

inline void cerror(const char *text, const char *file, const int line){
	fprintf(stderr,"Error %s\n", text);
	exit(1);
}

inline void cudasafe(int error, const char *message, const char *file, const int line){
	if (error != cudaSuccess) {
    cudaError cuerr = cudaGetLastError();
  	fprintf(stderr, "Error: %s in %s line %d\n", message, file, line); 
  	fprintf(stderr, "CUDA Error %d: %s\n", error, cudaGetErrorString(cuerr)); 
    exit(1);
  }
}

extern "C" {
	void gpu_info() {

  int nDevices;
	int i,j;
  size_t freemem;
  size_t total;

  cudasafe(cudaGetDeviceCount(&nDevices), "GetDeviceCount", __FILE__, __LINE__); 
  if (nDevices==0){
  	printf("No CUDA devices found"); 
		exit(-1);
	}
  printf("Number of CUDA devices: %d\n", nDevices); 

  for (i=0; i<nDevices; i++) {
    cudaDeviceProp prop;
    cudasafe(cudaGetDeviceProperties(&prop, i), "GetDeviceProperties", __FILE__, __LINE__);
    printf("*Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Compute Mode: %d\n", prop.computeMode);
    printf("  Can map host memory: %s\n", (prop.canMapHostMemory ? "Yes" : "No"));
    printf("  Major revision number: %d\n", prop.major);
    printf("  Minor revision number: %d\n", prop.minor);
    printf("  Total global memory: %u\n", prop.totalGlobalMem);
    printf("  Total shared memory per block: %u\n", prop.sharedMemPerBlock);
    printf("  Total registers per block: %d\n", prop.regsPerBlock);
    printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    printf("  Warp size: %d\n", prop.warpSize);
    printf("  Maximum memory pitch: %u\n", prop.memPitch);
    printf("  Maximum threads per block: %d\n", prop.maxThreadsPerBlock);
    for (j=0; j<3; ++j)
    	printf("  Maximum dimension %d of block: %d\n", j, prop.maxThreadsDim[j]);
    for (j=0; j<3; ++j)
    	printf("  Maximum dimension %d of grid: %d\n", j, prop.maxGridSize[j]);
    printf("  Clock rate (KHz): %d\n", prop.clockRate);
    printf("  Total constant memory: %u\n", prop.totalConstMem);
    printf("  Texture alignment: %u\n", prop.textureAlignment);
    printf("  Concurrent copy and execution: %s\n", (prop.deviceOverlap ? "Yes" : "No"));
    printf("  Number of multiprocessors: %d\n", prop.multiProcessorCount);
    printf("  Kernel execution timeout: %s\n\n", (prop.kernelExecTimeoutEnabled ? "Yes" : "No"));
    cudasafe(cudaSetDevice(i), "SetDevice", __FILE__, __LINE__);
    cudasafe(cudaMemGetInfo(&freemem,&total), "MemGetInfo", __FILE__, __LINE__);
  }

	}
}

