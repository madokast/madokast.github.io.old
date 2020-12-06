# -*- coding: utf-8 -*-

import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule

# 内核函数定义
mod = SourceModule("""
    // 需要从 CPU 端调用，让 GPU 执行的函数称为内核函数，需要用 __global__ 修饰
    __global__ void kernel_function(){
        printf("hello, world! -- from gpu\\n");
    }
""")

# 获取内核函数
kernel_function = mod.get_function("kernel_function")

# 调用内核函数
kernel_function(grid=(1,1,1), block=(1,1,1))
print("hello, pycuda! -- from python")

"""
注意：jupyter notebook 中没有打印 GPU 运行信息，可能是监听的 IO 流不对
在独立的 py 文件中，能收到 GPU 的打印
"""