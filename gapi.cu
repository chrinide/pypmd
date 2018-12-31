#include <stdio.h> 

#include "gapi.h"

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
