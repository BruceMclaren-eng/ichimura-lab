#include <stdio.h>

__global__ void vec_add(double *a, double *b, double *c, int n){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < n){
        c[i] = a[i] + b[i];
    }
}

int main(void){
    int n = 10;
    double a[10], b[10], c[10];
    double *d_a, *d_b, *d_c;

    // CPUの初期化
    for (int i = 0; i < n; i++){
        a[i] = 1.0;
        b[i] = 2.0;
    }

    // GPUメモリの確保
    cudaMalloc(&d_a, n * sizeof(double));
    cudaMalloc(&d_b, n * sizeof(double));
    cudaMalloc(&d_c, n * sizeof(double));

    // CPU→GPU転送
    cudaMemcpy(d_a, a, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, n * sizeof(double), cudaMemcpyHostToDevice);

    // カーネルの実行
    vec_add<<<1,n>>>(d_a, d_b, d_c, n);

    // GPU→CPU転送
    cudaMemcpy(c, d_c, n * sizeof(double), cudaMemcpyDeviceToHost);
    
    // 結果の表示
    for (int i = 0; i < n; i++){
        printf("c[%d] = %.1f\n", i, c[i]);
    }

    //GPUメモリの解放
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    return 0;
} 

