# try:
#     from books.cct.cctpy.cctpy import *
# except ModuleNotFoundError:
#     pass

# from cctpy import *

# 2020年12月8日 代码作废

GPU_ON: bool = True

try:
    import pycuda.autoinit
    import pycuda.driver as drv
    from pycuda.compiler import SourceModule
except ModuleNotFoundError as me:
    GPU_ON = False
    print("注意：没有安装 PYCUDA，无法使用 GPU 加速")


class GPU_ACCELERETE:
    COMPILED: bool = False  # 是否完成编译
    CUDA_MAGNETIC_FIELD_AT_CCT: Callable = None

    @staticmethod
    def __compile():
        if GPU_ACCELERETE.COMPILED:
            return  # 如果已经编译完毕，直接返回

        CUDA_GENERAL_CODE = """

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

    """

        CUDA_MAGNET_SOLO_CCT = """

// 计算单个 CCT 产生的磁场，winding 表示 CCT 路径离散点
// 因为 CUDA 只支持一维数组，winding[0]、winding[1]、winding[2]，表示第一个点
// winding[3]、winding[4]、winding[5] 表示第二个点
// length 表示点的数目
// 计算 CCT 在 p 点产生的磁场，返回值存入 ret 中
// 注意实际磁场还要乘上电流和 μ0/4π (=1e-7)
__global__ void magnet_solo_cct(float *winding, float *p, int *length, float *ret) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid == 0){
        vct_zero(ret);
        printf("\\%f,\\%f,\\%f\\n",p[0],p[1],p[2]);
    }

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
}

    """
        GPU_ACCELERETE.CUDA_MAGNETIC_FIELD_AT_CCT = SourceModule(
            CUDA_GENERAL_CODE + CUDA_MAGNET_SOLO_CCT
        ).get_function("magnet_solo_cct")

    @staticmethod
    def magnetic_field_at(magnet: Magnet, point: P3) -> P3:
        """
        magnet 在 point 处产生的磁场
        这个方法需要反复传输数据，速度比 CPU 慢
        """
        GPU_ACCELERETE.__compile()
        if isinstance(magnet, CCT):
            # point 转为局部坐标，并变成 numpy 向量
            p = magnet.local_coordinate_system.point_to_local_coordinate(
                point
            ).to_numpy_ndarry3_float32()
            length = int(magnet.dispersed_path3.shape[0])
            winding = magnet.dispersed_path3.flatten().astype(numpy.float32)
            ret = numpy.zeros((3,), dtype=numpy.float32)
            GPU_ACCELERETE.CUDA_MAGNETIC_FIELD_AT_CCT(
                drv.In(winding),
                drv.In(p),
                drv.In(numpy.array([length]).astype(numpy.int32)),
                drv.Out(ret),
                block=(512, 1, 1),
                grid=(256, 1),
            )
            print(p)
            return P3.from_numpy_ndarry(ret * magnet.current * 1e-7)