#include <stdio.h>

__global__ void kernel(double *x){
    x[threadIdx.x] = 2;
}

int main(){

	double *x;
  cudaMallocManaged(&x, sizeof(double)*2);
  cudaError_t error = cudaGetLastError();
  printf("%s\n", error);
  x[0] = 0;
  x[1] = 0;

  kernel<<<1, 2>>>(x);
  cudaDeviceSynchronize();

  printf("result = %f\n", x[0]);
  printf("result = %f\n", x[1]);

  cudaFree(x);
  return 0;
}
