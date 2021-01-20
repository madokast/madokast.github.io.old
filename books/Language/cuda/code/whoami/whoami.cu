#include <stdio.h>

 __global__ void report(){
    int i = blockIdx.x;
    int j = threadIdx.x;

    printf("My group id is %d, and my thread id is %d\n",i,j);
}

__global__ void report_in_detail(){
    int ix = blockIdx.x;
    int iy = blockIdx.y;
    int iz = blockIdx.z;

    int jx = threadIdx.x;
    int jy = threadIdx.y;
    int jz = threadIdx.z;

    printf("My group id is (%d,%d,%d), and my thread id is (%d,%d,%d)\n",ix,iy,iz,jx,jy,jz);
}

int main(){
    report_in_detail<<<3,2>>>();
    cudaThreadSynchronize(); // 同步标识。让 CPU 等待 GPU 运行结束
    printf("-----------\n");
    report_in_detail<<<dim3(1,1,3),dim3(1,2,1)>>>();
    return 0;
}