# -*- coding: utf-8 -*-

"""
直接分配显存
"""

import pycuda.driver as drv
import pycuda.autoinit
from pycuda.compiler import SourceModule

mod = SourceModule("""
    #include <stdio.h>
    __global__ void hello_kernel(){
        int N = 1024;
        int i;
        int *p = (int*)malloc(N*sizeof(int));
        for(i=0;i<N;i++){
            p[i] = i;
        }
        printf("%d",p[threadIdx.x]);
        free(p);
    }
""")

function = mod.get_function("hello_kernel")
function(block=(20, 1, 1),grid=(24,1))
print("hello, pycuda! -- from python in host")