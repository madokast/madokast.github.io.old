# -*- coding: utf-8 -*-

"""
CCT
"""

import pycuda.driver as drv
import pycuda.autoinit
from pycuda.compiler import SourceModule

mod = SourceModule("""
    #include <stdio.h>

    __global__ void cct_core(){
        printf("Hello,world!");
    }
""")

function = mod.get_function("cct_core")
function(block=(1, 1, 1),grid=(1,1))
print("hello, pycuda! -- from python in host")