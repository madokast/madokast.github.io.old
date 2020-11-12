# -*- coding: utf-8 -*-

import pycuda.driver as drv
import pycuda.autoinit
from pycuda.compiler import SourceModule

mod = SourceModule("""
    #include <stdio.h>
    __global__ void hello_kernel(){
        printf("hello, pycuda! -- from kernel in device");
    }
""")

function = mod.get_function("hello_kernel")
function(block=(1, 1, 1))
print("hello, pycuda! -- from python in host")
