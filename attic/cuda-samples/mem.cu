#include <stdio.h>


#define imin(a,b) (a<b?a:b)

const int N = 33 * 400 * 400;
const int threadsPerBlock = 1024;
const int blocksPerGrid =
            imin( 32, (N+threadsPerBlock-1) / threadsPerBlock );


__global__ void dot( int size, double *a, double *b, double *c ) {
    __shared__ double cache[threadsPerBlock];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;

    double   temp = 0;
    while (tid < size) {
        temp += a[tid] * b[tid];
        tid += blockDim.x * gridDim.x;
    }

    // set the cache values
    cache[cacheIndex] = temp;

    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = blockDim.x/2;
    while (i != 0) {
        if (cacheIndex < i)
            cache[cacheIndex] += cache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        c[blockIdx.x] = cache[0];
}


double malloc_test( int size ) {
    cudaEvent_t     start, stop;
    double           *a, *b, c, *partial_c;
    double           *dev_a, *dev_b, *dev_partial_c;
    float           elapsedTime;

    cudaEventCreate( &start );
    cudaEventCreate( &stop ) ;

    // allocate memory on the CPU side
    a = (double*)malloc( size*sizeof(double) );
    b = (double*)malloc( size*sizeof(double) );
    partial_c = (double*)malloc( blocksPerGrid*sizeof(double) );

    // allocate the memory on the GPU
    cudaMalloc( (void**)&dev_a,
                size*sizeof(double) ) ;
    cudaMalloc( (void**)&dev_b,
                size*sizeof(double) ) ;
    cudaMalloc( (void**)&dev_partial_c,
                              blocksPerGrid*sizeof(double) );

    // fill in the host memory with data
    for (int i=0; i<size; i++) {
        a[i] = i;
        b[i] = i*2;
    }

    cudaEventRecord( start, 0 ) ;
    // copy the arrays 'a' and 'b' to the GPU
    cudaMemcpy( dev_a, a, size*sizeof(double),
                cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size*sizeof(double),
                              cudaMemcpyHostToDevice ) ;

    dot<<<blocksPerGrid,threadsPerBlock>>>( size, dev_a, dev_b,
                                            dev_partial_c );
    // copy the array 'c' back from the GPU to the CPU
     cudaMemcpy( partial_c, dev_partial_c,
                 blocksPerGrid*sizeof(double),
                 cudaMemcpyDeviceToHost ) ;

     cudaEventRecord( stop, 0 ) ;
     cudaEventSynchronize( stop ) ;
     cudaEventElapsedTime( &elapsedTime,
                                        start, stop ) ;

    // finish up on the CPU side
    c = 0;
    for (int i=0; i<blocksPerGrid; i++) {
        c += partial_c[i];
    }

     cudaFree( dev_a ) ;
     cudaFree( dev_b ) ;
     cudaFree( dev_partial_c ) ;

    // free memory on the CPU side
    free( a );
    free( b );
    free( partial_c );

    // free events
     cudaEventDestroy( start  );
     cudaEventDestroy( stop  );

    printf( "Value calculated:  %f\n", c );

    return elapsedTime;
}


double cuda_host_alloc_test( int size ) {
    cudaEvent_t     start, stop;
    double           *a, *b, c, *partial_c;
    double           *dev_a, *dev_b, *dev_partial_c;
    float           elapsedTime;

    cudaEventCreate( &start ) ;
    cudaEventCreate( &stop ) ;

    // allocate the memory on the CPU
     cudaHostAlloc( (void**)&a,
                 size*sizeof(double),
                 cudaHostAllocWriteCombined |
                        cudaHostAllocMapped ) ;
     cudaHostAlloc( (void**)&b,
                 size*sizeof(double),
                 cudaHostAllocWriteCombined |
                        cudaHostAllocMapped ) ;
     cudaHostAlloc( (void**)&partial_c,
                 blocksPerGrid*sizeof(double),
                 cudaHostAllocMapped ) ;

    // find out the GPU pointers
     cudaHostGetDevicePointer( &dev_a, a, 0  );
     cudaHostGetDevicePointer( &dev_b, b, 0  );
     cudaHostGetDevicePointer( &dev_partial_c,
                                            partial_c, 0  );

    // fill in the host memory with data
    for (int i=0; i<size; i++) {
        a[i] = i;
        b[i] = i*2;
    }

     cudaEventRecord( start, 0 ) ;

    dot<<<blocksPerGrid,threadsPerBlock>>>( size, dev_a, dev_b,
                                            dev_partial_c );

  	cudaDeviceSynchronize();
    cudaEventRecord( stop, 0 ) ;
    cudaEventSynchronize( stop ) ;
    cudaEventElapsedTime( &elapsedTime,
                                        start, stop ) ;

    // finish up on the CPU side
    c = 0;
    for (int i=0; i<blocksPerGrid; i++) {
        c += partial_c[i];
    }

     cudaFreeHost( a ) ;
     cudaFreeHost( b ) ;
     cudaFreeHost( partial_c ) ;

    // free events
     cudaEventDestroy( start ) ;
     cudaEventDestroy( stop ) ;

    printf( "Value calculated:  %f\n", c );

    return elapsedTime;
}



// Main cuda function

int main() {
    printf("Hello from cuda \n");
    cudaDeviceProp  prop;
    int whichDevice;
     cudaGetDevice( &whichDevice ) ;
     cudaGetDeviceProperties( &prop, whichDevice ) ;
    if (prop.canMapHostMemory != 1) {
        printf( "Device can not map memory.\n" );
        return 0;
    }

    double           elapsedTime;

    cudaSetDeviceFlags( cudaDeviceMapHost );

    // try it with malloc
    elapsedTime = malloc_test( N );
    printf( "Time using cudaMalloc:  %3.1f ms\n",
            elapsedTime );

    // now try it with cudaHostAlloc
    elapsedTime = cuda_host_alloc_test( N );
    printf( "Time using cudaHostAlloc:  %3.1f ms\n",
            elapsedTime );
    return 0;
}
