#include <stdio.h>

__global__ void hello_kernel(void){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    printf("hello from thread %d\n",i);
}

int main(void){
    hello_kernel<<<2,5>>>();
    cudaDeviceSynchronize();
    return 0;
}