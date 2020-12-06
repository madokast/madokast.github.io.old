# -*- coding: utf-8 -*-

import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule
import numpy as np
import time

# 内核函数定义
mod = SourceModule("""
    // 需要从 CPU 端调用，让 GPU 执行的函数称为内核函数，需要用 __global__ 修饰
    __global__ void exhausted_word(float *a, int* length,float *ret){
        float sum = 0.0f;
        int len = *length;
        for(int i = 0;i < len;i++){
            for(int j = 0;j < 1000*120 + ((int)(sum-((int)sum)));j++){
                sum += sqrt(a[i]) + a[i]*a[i] + sqrt(a[i]*a[i]*a[i]);
            }
            sum /= 1000.0f;
        }
        *ret = sum;
    }
""")


exhausted_word = mod.get_function("exhausted_word")

a = np.linspace(0,100,101).astype(np.float32)
length = np.array([a.shape[0]],dtype=np.int32)
ret = np.empty(1,dtype=np.float32)

s = time.time()
exhausted_word(drv.In(a),drv.In(length),drv.Out(ret),grid=(1,1,1), block=(1,1,1))
print(f"CUDA ret={ret}, time={time.time()-s:.2f}秒")

s = time.time()
sum = 0.0
for i in range(a.shape[0]):
    for j in range(1000*120+int(sum-int(sum))):
        sum += np.sqrt(a[i]) + a[i]**2 + np.sqrt(a[i]**3);
    sum /= 1000.0
print(f"CPU ret={sum}, time={time.time()-s:.2f}秒")

