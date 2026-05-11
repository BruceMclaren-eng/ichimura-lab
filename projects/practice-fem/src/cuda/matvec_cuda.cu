#include <cuda_runtime.h>

// SpMV kernel：スレッドiが行iを計算
// データはすでにGPU上にある前提（転送なし）
__global__ void matvec_kernel(
    int n,
    const double *values,
    const int *col_idx,
    const int *row_ptr,
    const double *x,
    double *y)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    double tmp = 0.0;
    for (int k = row_ptr[i]; k < row_ptr[i+1]; k++) {
        tmp += values[k] * x[col_idx[k]];
    }
    y[i] = tmp;
}

// Fortranから呼ばれるラッパー
// GPU上のポインタをそのまま受け取る（転送なし）
extern "C" void matvec_cuda_device_(
    int *n,
    double *d_values,
    int *d_col_idx,
    int *d_row_ptr,
    double *d_x,
    double *d_y)
{
    int N = *n;
    int threads = 256;
    int blocks  = (N + threads - 1) / threads;
    matvec_kernel<<<blocks, threads>>>(N, d_values, d_col_idx, d_row_ptr, d_x, d_y);
    cudaDeviceSynchronize();
}