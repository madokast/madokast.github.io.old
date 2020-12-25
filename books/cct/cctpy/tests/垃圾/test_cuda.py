# -*- coding: utf-8 -*-

from cctpy import (
    BaseUtils,
    P2,
    P3,
    StraightLine2,
    Trajectory,
    Plot2,
    Plot3,
    CCT,
    LocalCoordinateSystem,
    MM,
)

import pycuda.autoinit
import pycuda.driver as drv
import numpy as np
from pycuda.compiler import SourceModule

mod = SourceModule(
    """
#include <stdio.h>
#include <math.h>
#include "cuda.h"

#define MM (0.001f)
#define DIM (3)
#define PI (3.1415927f)
#define X (0)
#define Y (1)
#define Z (2)


__device__ __forceinline__ void vct_cross(float *a, float *b, float *ret) {
    ret[X] = a[Y] * b[Z] - a[Z] * b[Y];
    ret[Y] = -a[X] * b[Z] + a[Z] * b[X];
    ret[Z] = a[X] * b[Y] - a[Y] * b[X];
}

__device__ __forceinline__ void vct_add_local(float *a_local, float *b) {
    a_local[X] += b[X];
    a_local[Y] += b[Y];
    a_local[Z] += b[Z];
}

__device__ __forceinline__ void vct_add(float *a, float *b, float *ret) {
    ret[X] = a[X] + b[X];
    ret[Y] = a[Y] + b[Y];
    ret[Z] = a[Z] + b[Z];
}

__device__ __forceinline__ void vct_dot_a_v(float a, float *v) {
    v[X] *= a;
    v[Y] *= a;
    v[Z] *= a;
}

__device__ __forceinline__ void vct_dot_a_v_ret(float a, float *v, float *ret) {
    ret[X] = v[X] * a;
    ret[Y] = v[Y] * a;
    ret[Z] = v[Z] * a;
}

__device__ __forceinline__ void vct_copy(float *src, float *des) {
    des[X] = src[X];
    des[Y] = src[Y];
    des[Z] = src[Z];
}

__device__ __forceinline__ float vct_len(float *v) {
    return sqrtf(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
}

__device__ __forceinline__ void vct_zero(float *v) {
    v[X] = 0.0f;
    v[Y] = 0.0f;
    v[Z] = 0.0f;
}

__device__ __forceinline__ void vct_sub(float *a, float *b, float *ret) {
    ret[X] = a[X] - b[X];
    ret[Y] = a[Y] - b[Y];
    ret[Z] = a[Z] - b[Z];
}

// 磁场计算 注意，这里计算的不是电流元的磁场，还需要乘以 电流 和 μ0/4π (=1e-7)
// p0 和 p1 构成电流元左右端点，计算电流元在 p 点产生的磁场，返回值写入 ret 中
__device__ void dB(float *p0, float *p1, float *p, float *ret) {
    float p01[DIM];
    float r[DIM];
    float rr;

    vct_sub(p1, p0, p01); // p01 = p1 - p0

    vct_add(p0, p1, r); // r = p0 + p1

    vct_dot_a_v(0.5f, r); // r = (p0 + p1)/2

    vct_sub(p, r, r); // r = p - r

    rr = vct_len(r); // rr = len(r)

    vct_cross(p01, r, ret); // ret = p01 x r

    rr = 1.0f / rr / rr / rr; // changed

    vct_dot_a_v(rr, ret); // rr . (p01 x r)
}

// 计算单个 CCT 产生的磁场，winding 表示 CCT 路径离散点
// 因为 CUDA 只支持一维数组，winding[0]、winding[1]、winding[2]，表示第一个点
// winding[3]、winding[4]、winding[5] 表示第二个点
// length 表示点的数目
// 计算 CCT 在 p 点产生的磁场，返回值存入 ret 中
// 注意实际磁场还要乘上电流和 μ0/4π (=1e-7)
__global__ void magnet_solo_cct(float *winding, float *p, int *length, float *ret) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid == 0)vct_zero(ret);

    __syncthreads();

    if (tid < *length - 1) {
        float *p0 = winding + tid * DIM;
        float *p1 = winding + (tid + 1) * DIM;
        float db[3];

        dB(p0, p1, p, db);

        atomicAdd(&ret[X], db[X]);
        atomicAdd(&ret[Y], db[Y]);
        atomicAdd(&ret[Z], db[Z]);
    }
}"""
)

magnet = mod.get_function("magnet")

cct = CCT(
    LocalCoordinateSystem.global_coordinate_system(),
    0.95,
    83 * MM + 15 * MM * 2,
    67.5,
    [30.0, 80.0, 90.0, 90.0],
    128,
    -9664,
    P2(0, 0),
    P2(128 * np.pi * 2, 67.5 / 180.0 * np.pi),
)

length = int(cct.dispersed_path3.shape[0])

print(f"len={length}")

winding = cct.dispersed_path3.flatten().astype(np.float32)


ret = np.empty((3,), dtype=np.float32)

p = np.array([0.0, 0.0, 0.0]).astype(np.float32)
magnet(
    drv.In(winding),
    drv.In(p),
    drv.In(np.array([length]).astype(np.int32)),
    drv.Out(ret),
    block=(512, 1, 1),
    grid=(250, 1),
)

print(ret)

print(ret * cct.current * 1e-7)

###################### time ###############
print("--------")
m = cct.magnetic_field_at_cpu(P3())
print(m)
m = cct.magnetic_field_at_gpu(P3())
print(m)

print("--------------------------------")

###################### time ###############
import time

times = 2
#################
s = time.time()
for x in np.linspace(0, 0.01, times):
    p = np.array([x, 0.0, 0.0]).astype(np.float32)
    magnet(
        drv.In(winding),
        drv.In(p),
        drv.In(np.array([length]).astype(np.int32)),
        drv.Out(ret),
        block=(512, 1, 1),
        grid=(250, 1),
    )
    print(P3.from_numpy_ndarry3(ret * cct.current * 1e-7), p)
print(f"CUDA_ONE d={time.time()-s}")

#################
s = time.time()
for x in np.linspace(0, 0.01, times):
    p = P3(x, 0, 0)
    m = cct.magnetic_field_at_cpu(p)
    print(m, p)
print(f"CPU-d={time.time()-s}")

#################
s = time.time()
# 主机数据
p_h = np.array([0.0, 0.0, 0.0]).astype(np.float32)
winding_h = winding  # winding = cct.dispersed_path3.flatten().astype(np.float32)
length_h = np.array([length]).astype(np.int32)
ret_h = np.empty((3,), dtype=np.float32)

# 设备数据
p_d = drv.mem_alloc(p_h.nbytes)
winding_d = drv.mem_alloc(winding_h.nbytes)
length_d = drv.mem_alloc(length_h.nbytes)
ret_d = drv.mem_alloc(ret_h.nbytes)

# 固定数据复制
drv.memcpy_htod(winding_d, winding_h)
drv.memcpy_htod(length_d, length_h)


for x in np.linspace(0, 0.01, times):
    # p 数据
    p_h[0] = x
    # 复制
    drv.memcpy_htod(p_d, p_h)

    magnet(
        winding_d,
        p_d,
        length_d,
        ret_d,
        block=(512, 1, 1),
        grid=(250, 1),
    )

    # 复制回来
    drv.memcpy_dtoh(ret_h, ret_d)
    print(P3.from_numpy_ndarry3(ret_h * cct.current * 1e-7), P3.from_numpy_ndarry3(p_h))
print(f"CUDA-TWO-d={time.time()-s}")