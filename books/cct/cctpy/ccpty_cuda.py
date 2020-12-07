# -*- coding: utf-8 -*-

"""
GPU CUDA 加速 cctpy 束流跟踪
"""

from typing import Literal
import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule
import numpy
import time

from scipy.optimize.optimize import main
from cctpy import *

try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass


class GPU_ACCELERATOR:
    FLOAT32: str = "FLOAT32"
    FLOAT64: str = "FLOAT64"

    def __init__(self, float_number_type: str = FLOAT32, block_dim_x: int = 1024) -> None:
        self.float_number_type = float_number_type

        if block_dim_x > 1024 or block_dim_x < 0:
            raise ValueError(
                f"block_dim_x 应 >=1 and <=1024 内取，不能是{block_dim_x}")
        if block_dim_x & (block_dim_x-1) != 0:
            raise ValueError(f"block_dim_x 应该取 2 的幂次，不能为{block_dim_x}")
        self.block_dim_x: int = int(block_dim_x)

        cuda_code_00_include = """
        #include <stdio.h>

        """

        cuda_code_01_float_type_define: str = None
        if float_number_type == GPU_ACCELERATOR.FLOAT32:
            cuda_code_01_float_type_define = """
            #define FLOAT32

            """
            self.numpy_dtype = numpy.float32
        elif float_number_type == GPU_ACCELERATOR.FLOAT64:
            cuda_code_01_float_type_define = """
            #define FLOAT64

            """
            self.numpy_dtype = numpy.float64
        else:
            raise ValueError(
                "float_number_type 必须是 GPU_ACCELERATOR.FLOAT32 或 GPU_ACCELERATOR.FLOAT64")

        # 头信息
        cuda_code_02_define = """
        #ifdef FLOAT32
        #define FLOAT float
        #else
        #define FLOAT double
        #endif


        #define DIM (3)
        #define X (0)
        #define Y (1)
        #define Z (2)
        #define BLOCK_DIM_X ({block_dim_x})
        """.format(block_dim_x=self.block_dim_x)

        # 向量运算内联函数
        cuda_code_03_vct_functions = """
        // 向量叉乘
        __device__ __forceinline__ void vct_cross(FLOAT *a, FLOAT *b, FLOAT *ret) {
            ret[X] = a[Y] * b[Z] - a[Z] * b[Y];
            ret[Y] = -a[X] * b[Z] + a[Z] * b[X];
            ret[Z] = a[X] * b[Y] - a[Y] * b[X];
        }

        // 向量原地加法
        __device__ __forceinline__ void vct_add_local(FLOAT *a_local, FLOAT *b) {
            a_local[X] += b[X];
            a_local[Y] += b[Y];
            a_local[Z] += b[Z];
        }

        // 向量原地加法
        __device__ __forceinline__ void vct6_add_local(FLOAT *a_local, FLOAT *b) {
            a_local[X] += b[X];
            a_local[Y] += b[Y];
            a_local[Z] += b[Z];
            a_local[X+DIM] += b[X+DIM];
            a_local[Y+DIM] += b[Y+DIM];
            a_local[Z+DIM] += b[Z+DIM];
        }

        // 向量加法
        __device__ __forceinline__ void vct_add(FLOAT *a, FLOAT *b, FLOAT *ret) {
            ret[X] = a[X] + b[X];
            ret[Y] = a[Y] + b[Y];
            ret[Z] = a[Z] + b[Z];
        }

        // 向量加法
        __device__ __forceinline__ void vct6_add(FLOAT *a, FLOAT *b, FLOAT *ret) {
            ret[X] = a[X] + b[X];
            ret[Y] = a[Y] + b[Y];
            ret[Z] = a[Z] + b[Z];
            ret[X+DIM] = a[X+DIM] + b[X+DIM];
            ret[Y+DIM] = a[Y+DIM] + b[Y+DIM];
            ret[Z+DIM] = a[Z+DIM] + b[Z+DIM];
        }

        // 向量*常数，原地操作
        __device__ __forceinline__ void vct_dot_a_v(FLOAT a, FLOAT *v) {
            v[X] *= a;
            v[Y] *= a;
            v[Z] *= a;
        }

        // 向量*常数，原地操作
        __device__ __forceinline__ void vct6_dot_a_v(FLOAT a, FLOAT *v) {
            v[X] *= a;
            v[Y] *= a;
            v[Z] *= a;
            v[X+DIM] *= a;
            v[Y+DIM] *= a;
            v[Z+DIM] *= a;
        }

        // 向量*常数
        __device__ __forceinline__ void vct_dot_a_v_ret(FLOAT a, FLOAT *v, FLOAT *ret) {
            ret[X] = v[X] * a;
            ret[Y] = v[Y] * a;
            ret[Z] = v[Z] * a;
        }

        // 向量*常数
        __device__ __forceinline__ void vct6_dot_a_v_ret(FLOAT a, FLOAT *v, FLOAT *ret) {
            ret[X] = v[X] * a;
            ret[Y] = v[Y] * a;
            ret[Z] = v[Z] * a;
            ret[X+DIM] = v[X+DIM] * a;
            ret[Y+DIM] = v[Y+DIM] * a;
            ret[Z+DIM] = v[Z+DIM] * a;
        }

        __device__ __forceinline__ FLOAT vct_dot_v_v(FLOAT *v,FLOAT *w){
            return v[X] * w[X] + v[Y] * w[Y] + v[Z] * w[Z];
        }

        // 向量拷贝赋值
        __device__ __forceinline__ void vct_copy(FLOAT *src, FLOAT *des) {
            des[X] = src[X];
            des[Y] = src[Y];
            des[Z] = src[Z];
        }

        // 求向量长度
        __device__ __forceinline__ FLOAT vct_len(FLOAT *v) {

            #ifdef FLOAT32
            return sqrtf(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
            #else
            return sqrt(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
            #endif
        }

        // 将矢量 v 置为 0
        __device__ __forceinline__ void vct_zero(FLOAT *v) {
            v[X] = 0.0;
            v[Y] = 0.0;
            v[Z] = 0.0;
        }

        // 打印矢量，一般用于 debug
        __device__ __forceinline__ void vct_print(FLOAT *v) {
            #ifdef FLOAT32
            printf("%.15f, %.15f, %.15f\\n", v[X], v[Y], v[Z]);
            #else
            printf("%.15lf, %.15lf, %.15lf\\n", v[X], v[Y], v[Z]);
            #endif
        }

        // 矢量减法
        __device__ __forceinline__ void vct_sub(FLOAT *a, FLOAT *b, FLOAT *ret) {
            ret[X] = a[X] - b[X];
            ret[Y] = a[Y] - b[Y];
            ret[Z] = a[Z] - b[Z];
        }

        """

        cuda_code_04_dB = """
        // 计算电流元在 p 点产生的磁场
        // 其中 p0 表示电流元的位置
        // kl 含义见下
        // 返回值放在 ret 中
        // 
        // 原本电流元的计算公式如下：
        // dB = (miu0/4pi) * Idl × r / (r^3)
        // 其中 r = p - p0，p0 是电流元的位置
        // 
        // 如果考虑极小一段电流（起点s0，终点s1）则产生的磁场为
        // ΔB = (miu0/4pi) * I * (s1-s2)*r / (r^3)
        // 同样的，r = p - p0，p0 = (s1+s2)/2
        //
        // 因为 (miu0/4pi) * I * (s1-s2) 整体已知，所以提前计算为 kl
        // p0 提前已知，即 (s1+s2)/2，也提前给出
        // 这样可以减少无意义的重复计算
        //
        // 补充：坐标均是全局坐标
        __device__ __forceinline__ void dB(FLOAT *kl, FLOAT *p0, FLOAT *p, FLOAT *ret){
            FLOAT r[DIM];
            FLOAT rr;

            vct_sub(p, p0, r); // r = p - p0

            rr = vct_len(r); // rr = abs(r)

            rr = rr*rr*rr; // rr = rr^3

            vct_cross(kl, r, ret); // ret = kl × r

            vct_dot_a_v(1.0/rr, ret); // ret = (kl × r)/(rr^3)
        }

        // 计算所有的电流元在 p 点产生的磁场
        // number 表示电流元数目
        // kls 每 DIM = 3 组表示一个 kl
        // p0s 每 DIM = 3 组表示一个 p0
        __device__ void current_element_B(FLOAT *kls, FLOAT *p0s, int number, FLOAT *p, FLOAT *ret){
            int tid = threadIdx.x; // 0-1023 (decide by BLOCK_DIM_X)
            FLOAT db[DIM];
            __shared__ FLOAT s_dbs[DIM*BLOCK_DIM_X];

            vct_zero(s_dbs + tid*DIM);

            // 计算每个电流元产生的磁场
            for(int i = tid*DIM; i < number*DIM; i += BLOCK_DIM_X*DIM){
                dB(kls + i, p0s + i, p, db);
                vct_add_local(s_dbs + tid*DIM, db);
            }
            
            // 规约求和（from https://www.bilibili.com/video/BV15E411x7yT）
            for(int step = BLOCK_DIM_X>>1; step >= 1; step>>=1){
                __syncthreads(); // 求和前同步
                if(tid<step) vct_add_local(s_dbs + tid * DIM, s_dbs + (tid + step) * DIM);
            }

            if(tid == 0) vct_copy(s_dbs, ret);
        }

        """

        cuda_code_05_QS = """
        // 计算 QS 在 p 点产生的磁场
        // origin xi yi zi 分别是 QS 的局部坐标系
        // 这个函数只需要单线程计算
        __device__ __forceinline__ void magnet_at_qs(FLOAT *origin, FLOAT *xi, 
                FLOAT *yi, FLOAT *zi, FLOAT length, FLOAT gradient, FLOAT second_gradient, FLOAT aper_r, FLOAT *p, FLOAT* ret){
            FLOAT temp1[DIM];
            FLOAT temp2[DIM];

            vct_sub(p, origin, temp1); // temp1 = p - origin
            temp2[X] = vct_dot_v_v(xi, temp1);
            temp2[Y] = vct_dot_v_v(yi, temp1);
            temp2[Z] = vct_dot_v_v(zi, temp1); // 这时 temp2 就是全局坐标 p 点在 QS 局部坐标系中的坐标

            vct_zero(ret);

            if(temp2[Z]<0 || temp2[Z]>length){
                return; // 无磁场
            }else{
                if(
                    temp2[X] > aper_r ||
                    temp2[X] < -aper_r ||
                    temp2[Y] > aper_r ||
                    temp2[Y] < -aper_r ||
                    #ifdef FLOAT32
                    sqrtf(temp2[X]*temp2[X]+temp2[Y]*temp2[Y]) > aper_r
                    #else
                    sqrt(temp2[X]*temp2[X]+temp2[Y]*temp2[Y]) > aper_r
                    #endif
                ){
                    return; // 无磁场
                }else{
                    temp1[X] = gradient * temp2[Y] + second_gradient * (temp2[X] * temp2[Y]);
                    temp1[Y] = gradient * temp2[X] + 0.5 * second_gradient * (temp2[X] * temp2[X] - temp2[Y] * temp2[Y]);

                    vct_dot_a_v_ret(temp1[X], xi, ret);
                    vct_dot_a_v_ret(temp1[Y], yi, temp2);
                    vct_add_local(ret, temp2);
                }
            }
        }
        
        """

        cuda_code_06_runge_kutta4 = """
        // runge_kutta4 代码和 cctpy 中的 runge_kutta4 一模一样
        // Y0 数组长度为 6
        // Y0 会发生变化，既是输入也是输出
        // 为了分析包络等，会出一个记录全部 YO 的函数
        // 这个函数也是单线程运行
        __device__ void runge_kutta4(FLOAT t0, FLOAT t_end, FLOAT *Y0, void (*call)(FLOAT,FLOAT*,FLOAT*), FLOAT dt){
            #ifdef FLOAT32
            int number = (int)(ceilf((t_end - t0) / dt));
            #else
            int number = (int)(ceil((t_end - t0) / dt));
            #endif

            dt = (t_end - t0) / ((FLOAT)(number));
            FLOAT k1[DIM*2];
            FLOAT k2[DIM*2];
            FLOAT k3[DIM*2];
            FLOAT k4[DIM*2];
            FLOAT temp[DIM*2];

            for(int ignore = 0; ignore < number; ignore++){
                (*call)(t0, Y0, k1);

                vct6_dot_a_v_ret(dt / 2., k1, temp); // temp = dt / 2 * k1
                vct6_add_local(temp, Y0); // temp =  Y0 + temp
                (*call)(t0 + dt / 2., temp, k2);


                vct6_dot_a_v_ret(dt / 2., k2, temp); // temp = dt / 2 * k2
                vct6_add_local(temp, Y0); // temp =  Y0 + temp
                (*call)(t0 + dt / 2., temp, k3);

                vct6_dot_a_v_ret(dt, k3, temp); // temp = dt * k3
                vct6_add_local(temp, Y0); // temp =  Y0 + temp
                (*call)(t0 + dt, temp, k4);

                t0 += dt;
                
                vct6_add(k1, k4, temp); // temp = k1 + k4
                vct6_dot_a_v(2.0, k2);
                vct6_dot_a_v(2.0, k3);
                vct6_add(k2, k3, k1); // k1 已经没用了，所以装 k1 = k2 + k3
                vct6_add_local(temp, k1);
                vct6_dot_a_v(dt / 6.0, temp);
                vct6_add_local(Y0, temp);
                // Y0 += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
            }

            
        }
        """

        cuda_code_07_run_only = """
        // 单个粒子跟踪
        // FLOAT *kls, FLOAT *p0s, int number 电流元相关
        // FLOAT *qs_data, int qs_number
        // 这个函数也是单线程运行
        __device__ void track_particle(Float* particle, Float length, Float footstep,
            FLOAT *kls, FLOAT *p0s, int current_element_number, FLOAT *qs_data, int qs_number){
            // TODO 
            // 确定参数

        }
        
        """


        self.cuda_code: str = (
            cuda_code_00_include +
            cuda_code_01_float_type_define +
            cuda_code_02_define +
            cuda_code_03_vct_functions +
            cuda_code_04_dB +
            cuda_code_05_QS + 
            cuda_code_06_runge_kutta4
        )

    def vct_length(self, p3: P3):
        """
        测试用函数，计算矢量长度
        示例：
        ga32 = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT32)
        ga64 = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT64)
        v = P3(1,1,1)
        print(f"diff={ga32.vct_length(v) - v.length()}") # diff=-3.1087248775207854e-08
        print(f"diff={ga64.vct_length(v) - v.length()}") # diff=0.0
        """
        code = """
        __global__ void vl(FLOAT* v, FLOAT* ret){
            *ret = vct_len(v);
        }

        """

        mod = SourceModule(self.cuda_code + code)

        vl = mod.get_function("vl")

        ret = numpy.empty((1,), dtype=self.numpy_dtype)

        vl(drv.In(p3.to_numpy_ndarry3(numpy_dtype=self.numpy_dtype)),
           drv.Out(ret), grid=(1, 1, 1), block=(1, 1, 1))

        return float(ret[0])

    def vct_print(self, p3: P3):
        """
        测试用函数，打印矢量
        示例：
        ga32 = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT32)
        ga64 = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT64)
        v = P3(1/3, 1/6, 1/7)
        ga32.vct_print(v)
        ga64.vct_print(v)
        >>>
        0.333333343267441, 0.166666671633720, 0.142857149243355
        0.333333333333333, 0.166666666666667, 0.142857142857143
        """
        code = """
        __global__ void vp(FLOAT* v){
            vct_print(v);
        }

        """

        mod = SourceModule(self.cuda_code + code)

        vp = mod.get_function("vp")

        vp(drv.In(p3.to_numpy_ndarry3(numpy_dtype=self.numpy_dtype)),
           grid=(1, 1, 1), block=(1, 1, 1))

    def current_element_B(self, kls: numpy.ndarray, p0s: numpy.ndarray, number: int, p: P3):
        """
        计算电流元集合，在 p 点产生的磁场
        对比代码如下：
        -------------------
        bl = HUST_SC_GANTRY.beamline
        cct:CCT = bl.magnets[15]
        p = bl.trajectory.point_at(HUST_SC_GANTRY.beamline_length_part1+HUST_SC_GANTRY.DL2+0.5).to_p3() + P3(1E-3,1E-4,1E-5);
        print(cct.magnetic_field_at(p))

        ga = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT64)

        kls,p0s = cct.global_current_elements_and_elementary_current_positions(numpy_dtype=numpy.float64)

        b = ga.current_element_B(
            kls.flatten(),
            p0s.flatten(),
            cct.total_disperse_number,
            p,
        )

        print(b)
        -------------------
        >>>
        float32下
        [-1.3026414984968506, -0.34195788415622125, 1.249137357827171]
        [-1.302641749382019,  -0.34195759892463684, 1.2491374015808105] 1.8368872515000503
        （小数 6 位相同）
        float64下
        [-1.3026414984968506, -0.34195788415622125, 1.249137357827171]
        [-1.302641498496853,  -0.3419578841562222,  1.24913735782717] 1.836887096928434
        （小数 14 位相同，看来绝不会出现电流元数目多一少一的错误）
        """

        code = """
        __global__ void ce(FLOAT *kls, FLOAT *p0s, int* number, FLOAT *p, FLOAT *ret){
            current_element_B(kls,p0s,*number,p,ret);
        }

        """

        mod = SourceModule(self.cuda_code + code)

        ce = mod.get_function("ce")

        ret = numpy.empty((3,), dtype=self.numpy_dtype)

        ce(drv.In(kls.astype(self.numpy_dtype)),
           drv.In(p0s.astype(self.numpy_dtype)),
           drv.In(numpy.array([number], dtype=numpy.int32)),
           drv.In(p.to_numpy_ndarry3(numpy_dtype=self.numpy_dtype)),
           drv.Out(ret),
           grid=(1, 1, 1), block=(self.block_dim_x, 1, 1))

        return P3.from_numpy_ndarry(ret)

    def magnet_at_qs(self, qs_data, p3: P3):
        """
        qs 磁铁在 p 点产生的磁场
        p 点是 全局坐标点
        测试：
        ---------------------------
            bl = HUST_SC_GANTRY.beamline
        qs: QS = bl.magnets[23]
        p = bl.trajectory.point_at(HUST_SC_GANTRY.beamline_length_part1+HUST_SC_GANTRY.DL2 +
                                1.19+HUST_SC_GANTRY.GAP1+HUST_SC_GANTRY.qs3_length/2).to_p3() + P3(10*MM, 10*MM, 10*MM)
        print(qs.magnetic_field_at(p))
        ga = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT32)

        qs_data = numpy.array(
            qs.local_coordinate_system.location.to_list() + qs.local_coordinate_system.XI.to_list()
            + qs.local_coordinate_system.YI.to_list() + qs.local_coordinate_system.ZI.to_list() 
            + [qs.length, qs.gradient, qs.second_gradient, qs.aperture_radius]
            ,dtype=numpy.float32)

        b =  ga.magnet_at_qs(
            qs_data=qs_data,
            p3=p
        )

        print(b)
        ---------------------------
        >>> float32
        [-0.03274739751857501, -0.07905921122176701, -0.09954070791478388]
        [-0.03274737671017647, -0.0790591612458229, -0.09954022616147995]
        >>> float64
        [-0.03274739751857501, -0.07905921122176701, -0.09954070791478388]
        [-0.03274739751857501, -0.07905921122176701, -0.09954070791478384]
        """

        code = """
        __global__ void mq(FLOAT *qs_data, FLOAT *p, FLOAT *ret){
            magnet_at_qs(
                qs_data, // origin
                qs_data + 3, //xi
                qs_data + 6, //yi
                qs_data + 9, //zi
                *(qs_data + 12), // len
                *(qs_data + 13), // g
                *(qs_data + 14), // sg
                *(qs_data + 15), // aper r
                p, ret
            );
        }

        """

        mod = SourceModule(self.cuda_code + code)

        mq = mod.get_function("mq")

        ret = numpy.empty((3,), dtype=self.numpy_dtype)

        mq(drv.In(qs_data.astype(self.numpy_dtype)),
        drv.In(p.to_numpy_ndarry3(numpy_dtype=self.numpy_dtype)),
        drv.Out(ret),
        grid=(1, 1, 1), block=(1, 1, 1)
        )

        return P3.from_numpy_ndarry(ret)


if __name__ == "__main__":
    bl = HUST_SC_GANTRY.beamline
    qs: QS = bl.magnets[23]
    p = bl.trajectory.point_at(HUST_SC_GANTRY.beamline_length_part1+HUST_SC_GANTRY.DL2 +
                               1.19+HUST_SC_GANTRY.GAP1+HUST_SC_GANTRY.qs3_length/2).to_p3() + P3(10*MM, 10*MM, 10*MM)
    print(qs.magnetic_field_at(p))
    ga = GPU_ACCELERATOR(float_number_type=GPU_ACCELERATOR.FLOAT64)

    qs_data = numpy.array(
        qs.local_coordinate_system.location.to_list() + qs.local_coordinate_system.XI.to_list()
        + qs.local_coordinate_system.YI.to_list() + qs.local_coordinate_system.ZI.to_list() 
        + [qs.length, qs.gradient, qs.second_gradient, qs.aperture_radius]
        ,dtype=numpy.float64)

    b =  ga.magnet_at_qs(
        qs_data=qs_data,
        p3=p
    )

    print(b)
