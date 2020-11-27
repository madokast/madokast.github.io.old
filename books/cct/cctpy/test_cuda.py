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

__global__ void magnet(float *winding, float *p, int *length, float *ret) {
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

###################### time ###############
import time

s = time.time()
for x in np.linspace(0,1,100):
    p = np.array([x, 0.0, 0.0]).astype(np.float32)
    magnet(
        drv.In(winding),
        drv.In(p),
        drv.In(np.array([length]).astype(np.int32)),
        drv.Out(ret),
        block=(512, 1, 1),
        grid=(250, 1),
    )
print(f"d={time.time()-s}")

s = time.time()
for x in np.linspace(0,1,100):
    m = cct.magnetic_field_at(P3(x,0,0))
print(f"d={time.time()-s}")