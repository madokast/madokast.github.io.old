#include <stdio.h>

 __global__ void task_allocate(){
    int i = threadIdx.x;
    if(i==0){ // 计算 1+2+3+...+100
        int n = 100;
        int sum = 0;
        while(n>0){
            sum+=n;
            n--;
        }
        printf("1+2+...+100 = %d\n",sum);
    }else if(i==1){ // 计算 10 的阶乘
        int n = 10;
        int factor = 1;
        while(n>0){
            factor*=n;
            n--;
        }
        printf("10! = %d\n",factor);
    }
}

int main(){
    task_allocate<<<1,2>>>();
    return 0;
}