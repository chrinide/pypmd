#define N 512
#define BLOCK_DIM 512
__global__ void matrixAdd (int *a, int *b, int *c);
int main() {
 int a[N][N], b[N][N], c[N][N];
 int *dev_a, *dev_b, *dev_c;
 int size = N * N * sizeof(int);
 // initialize a and b with real values (NOT SHOWN)
 cudaMalloc((void**)&dev_a, size);
 cudaMalloc((void**)&dev_b, size);
 cudaMalloc((void**)&dev_c, size);
 cudaMemcpy(dev_a, a, size, cudaMemcpyHostToDevice);
 cudaMemcpy(dev_b, b, size, cudaMemcpyHostToDevice);
 dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
 dim3 dimGrid((int)ceil(N/dimBlock.x),(int)ceil(N/dimBlock.y));
 matrixAdd<<<dimGrid,dimBlock>>>(dev_a,dev_b,dev_c);
 cudaMemcpy(c, dev_c, size, cudaMemcpyDeviceToHost);
 cudaFree(dev_a); cudaFree(dev_b); cudaFree(dev_c);
}
__global__ void matrixAdd (int *a, int *b, int *c) {
 int col = blockIdx.x * blockDim.x + threadIdx.x;
 int row = blockIdx.y * blockDim.y + threadIdx.y;
 int idx = col + row * N;
 if (col < N && row < N) {
 c[idx] = a[idx] + b[idx];
 }
}
__global__ void matrixMult (int *a, int *b, int *c, int width) {
 int k, sum = 0;
 int col = threadIdx.x + blockDim.x * blockIdx.x;
 int row = threadIdx.y + blockDim.y * blockIdx.y;
 if(col < width && row < width) {
 for (k = 0; k < width; k++)
 sum += a[row * width + k] * b[k * width + col];
 c[row * width + col] = sum;
 }
}
