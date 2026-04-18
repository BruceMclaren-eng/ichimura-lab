cat << 'EOF' > device_info.cu
#include <stdio.h>
int main(void) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("最大スレッド数/ブロック : %d\n", prop.maxThreadsPerBlock);
    printf("最大ブロック数(x)       : %d\n", prop.maxGridSize[0]);
    printf("SM数                    : %d\n", prop.multiProcessorCount);
    printf("SM当たり最大スレッド数  : %d\n", prop.maxThreadsPerMultiProcessor);
    printf("総最大スレッド数        : %d\n", prop.multiProcessorCount * prop.maxThreadsPerMultiProcessor);
    return 0;
}
EOF