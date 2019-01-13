#include <stdio.h> 

void cudasafe(int error, const char *message, const char *file, const int line){
	if (error != cudaSuccess) {
  	fprintf(stderr, "CUDA Error: %s : %i. In %s line %d\n", message, error, file, line); 
    exit(-1);
  }
}

void cuda_error(FILE *fp, const char *message){
	cudaError_t error = cudaGetLastError();
  if (error!=cudaSuccess) {
     fprintf(fp,"\n  ERROR: %s: %s\n\n", message, cudaGetErrorString(error) );
     fflush(fp);
     exit(-1);
  }
}

int main() {

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
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
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
    cudaSetDevice(i);
    cuda_error(stdout,"cudasetdevice");
    cudaMemGetInfo(&freemem,&total);
    cuda_error(stdout,"cudamemgetinfo");
  }

	return 0;
}
