# -*- coding: utf-8 -*-

from numpy.core.defchararray import add
from numpy.core.shape_base import block
import pycuda.autoinit
import pycuda.driver as drv
import numpy as np
from pycuda.compiler import SourceModule
import time

mod = SourceModule(
    """
    __global__ void timer(float *a){
        int max = 1000*1000;
        int i;
        float temp = 0.0f;
        for(i=0;i<max;i++){
            temp += + sinf(i);
        }
        a[0] = temp;
    }
"""
)


for i in range(5):
    timer = mod.get_function("timer")

    N = 512
    B = drv.Device(0).get_attributes()[pycuda._driver.device_attribute.MULTIPROCESSOR_COUNT]

    h_a = np.linspace(0.1, 0.5, 1).astype(np.float32)

    start = time.time()

    timer(drv.InOut(h_a), block=(N, 1, 1), grid=(B, 1))

    print(h_a[0])

    end = time.time()

    print(f"time={end-start}")

# temp = np.array([0.0],dtype=np.float64)
# for i in range(1000*1000):
#     temp += np.sin(i)
# print(temp)