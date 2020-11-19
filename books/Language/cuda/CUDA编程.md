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

## 安装

- 安装 visual studio 2019，安装时选择"使用C++的桌面开发"，之后在安装目录搜索64位的 C/C++ 编译器 cl.exe，把目录加入 path 环境变量。

```
典型 cl.exe 目录
C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\bin\Hostx64\x64
```

- 安装 NVIDIA 驱动程序，https://www.nvidia.cn/geforce/drivers/ ，完成驱动更新。

- 在 NVIDIA control panel 软件中，查看 CUDA 的版本，方法如下。图片中 CUDA 版本为11.1
  
<img src="./img/NVIDIA_control_panel_main.jpg"></img>

<img src="./img/NVIDIA_control_panel_sys_info.jpg"></img>

- 安装 CUDA 开发工具，https://developer.nvidia.com/cuda-zone ，选择对应的版本。历史版本见 https://developer.nvidia.com/cuda-toolkit-archive 

- 安装 pycuda，`pip install -i https://pypi.tuna.tsinghua.edu.cn/simple pycuda`

## GPU 架构

- 物理架构

GPU 硬件的核心是 SM，streaming multiprocessor 流式多处理器。SM 包含 CUDA 核心、共享内存、寄存器等，可以并发的执行数百个线程。当 kernel 被执行时，grid 中的线程被分配到一个 SM 上。

SM 采用 SIMT，single instruction multipthread，单指令多线程结构，基本的执行单位是线程束 wrap，线程束包含 32 个线程。线程束上的线程执行相同的指令，但是因为有各自的寄存器，所以执行相同指令，但是控制的地址和计算结果可以不同。

线程束中的线程遇到了分支结构时，不进入分支的线程只能等待。所以线程束中出现程序分支会导致性能下降。

- 编程架构 网格、线程块

一个 kernel 的所有线程，在物理层不一定是全并发的。因为线程束 wrap 是 32，所以 block 的大小一般设置为 32 的整数倍。

## 内置变量 blockidx threadidx 三维

## 规约算法 reduction

- 定义，将多个变量的通过某种运算变成一个值，如数组求和、向量内积。单线程一般使用循环。

- 并行下，可以使用 树型并行 进行规约。

<img src="./img/规约算法.jpg"></img>

- 同步操作，需要同步等待同一行的规约完成。但是 CUDA 不支持全局同步，因为当程序块数目大于 SM 数目时，会发生死锁，没有空闲的物理线程去执行剩余的同行的规约运算。

## 注意事项

不要多个线程访问同一个地址！这样会发生冲突，导致效率低下。