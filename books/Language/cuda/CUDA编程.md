# CUDA 编程

## 基本概念

- CUDA（compute unified device architecture 计算统一设备架构），英伟达 NVIDIA 开发的并行计算模型，支持 NVIDIA GPU卡（图形处理单元）

- GPU 最开始用于渲染图象，也可用于其他数学计算。

- 和 CPU 相比，GPU 的控制部分简单，数据计算部分复杂。

- CUDA 程序一般流程

```
1. 设备（GPU）分配内存
2. 数据从主机复制到设备
3， 指定并行度，启动内核计算
4. 数据复制回主机
5. 释放内存
```

- 开发 CUDA 程序，需要 GPU 编译器