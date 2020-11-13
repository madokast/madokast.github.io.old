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

- SM 流多处理器，下有若干 SP 核心。

- 每 32 个线程共享指令。

## 注意事项

不要多个线程访问同一个地址！这样会发生冲突。