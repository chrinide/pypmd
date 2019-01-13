#include <stdio.h>

__global__ void MyKernel(int *a, int *b, int *c, int N) 
{ 
    int idx = threadIdx.x + blockIdx.x * blockDim.x; 

    if (idx < N) { c[idx] = a[idx] + b[idx]; } 
} 

int main() 
{ 
  const int N = 1000000;

  int blockSize;      // The launch configurator returned block size 
  int minGridSize;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch 
  int gridSize;       // The actual grid size needed, based on input size 

  int *h_vec1 = (int*) malloc(N*sizeof(int));
  int *h_vec2 = (int*) malloc(N*sizeof(int));
  int *h_vec3 = (int*) malloc(N*sizeof(int));
  int *h_vec4 = (int*) malloc(N*sizeof(int));

  int *d_vec1; cudaMalloc((void**)&d_vec1, N*sizeof(int));
  int *d_vec2; cudaMalloc((void**)&d_vec2, N*sizeof(int));
  int *d_vec3; cudaMalloc((void**)&d_vec3, N*sizeof(int));

  for (int i=0; i<N; i++) {
  	h_vec1[i] = 10;
    h_vec2[i] = 20;
    h_vec4[i] = h_vec1[i] + h_vec2[i];
  }

  cudaMemcpy(d_vec1, h_vec1, N*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_vec2, h_vec2, N*sizeof(int), cudaMemcpyHostToDevice);

  float time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, MyKernel, 0, N); 

  // Round up according to array size 
  gridSize = (N + blockSize - 1) / blockSize; 

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("Occupancy calculator elapsed time:  %3.3f ms \n", time);

  cudaEventRecord(start, 0);

  MyKernel<<<gridSize, blockSize>>>(d_vec1, d_vec2, d_vec3, N); 

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("Kernel elapsed time:  %3.3f ms \n", time);

  printf("Blocksize %i\n", blockSize);

  cudaMemcpy(h_vec3, d_vec3, N*sizeof(int), cudaMemcpyDeviceToHost);

  for (int i=0; i<N; i++) {
  	if (h_vec3[i] != h_vec4[i]) { 
			printf("Error at i = %i! Host = %i; Device = %i\n", i, h_vec4[i], h_vec3[i]);
			return 0; 
		}
  }

  printf("Test passed\n");

	return 0;

}
