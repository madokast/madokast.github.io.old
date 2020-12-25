#include <stdio.h>


// 定义为 32 位浮点数模式
#define FLOAT32


// 根据定义的浮点数模式，将 FLOAT 宏替换为 float 或 double
#ifdef FLOAT32
#define FLOAT float
#else
#define FLOAT double
#endif

// 维度 三维
#define DIM (3)
// 维度索引 0 1 2 表示 X Y Z，这样对一个数组取值，看起来清晰一些
#define X (0)
#define Y (1)
#define Z (2)
// 粒子参数索引 (px0, py1, pz2, vx3, vy4, vz5, rm6 相对质量, e7 电荷量, speed8 速率, distance9 运动距离)
#define PARTICLE_DIM (10)
#define PX (0)
#define PY (1)
#define PZ (2)
#define VX (3)
#define VY (4)
#define VZ (5)
#define RM (6)
#define E (7)
#define SPEED (8)
#define DISTANCE (9)

// 块线程数目
#define BLOCK_DIM_X (1024)
#define QS_DATA_LENGTH (16)
#define MAX_CURRENT_ELEMENT_NUMBER (240000)

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
    a_local[X + DIM] += b[X + DIM];
    a_local[Y + DIM] += b[Y + DIM];
    a_local[Z + DIM] += b[Z + DIM];
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
    ret[X + DIM] = a[X + DIM] + b[X + DIM];
    ret[Y + DIM] = a[Y + DIM] + b[Y + DIM];
    ret[Z + DIM] = a[Z + DIM] + b[Z + DIM];
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
    v[X + DIM] *= a;
    v[Y + DIM] *= a;
    v[Z + DIM] *= a;
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
    ret[X + DIM] = v[X + DIM] * a;
    ret[Y + DIM] = v[Y + DIM] * a;
    ret[Z + DIM] = v[Z + DIM] * a;
}

__device__ __forceinline__ FLOAT vct_dot_v_v(FLOAT *v, FLOAT *w) {
    return v[X] * w[X] + v[Y] * w[Y] + v[Z] * w[Z];
}

// 向量拷贝赋值
__device__ __forceinline__ void vct_copy(FLOAT *src, FLOAT *des) {
    des[X] = src[X];
    des[Y] = src[Y];
    des[Z] = src[Z];
}

// 向量拷贝赋值
__device__ __forceinline__ void vct6_copy(FLOAT *src, FLOAT *des) {
    des[X] = src[X];
    des[Y] = src[Y];
    des[Z] = src[Z];
    des[X + DIM] = src[X + DIM];
    des[Y + DIM] = src[Y + DIM];
    des[Z + DIM] = src[Z + DIM];
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
    printf("%.15f, %.15f, %.15f\n", v[X], v[Y], v[Z]);
#else
    printf("%.15lf, %.15lf, %.15lf\n", v[X], v[Y], v[Z]);
#endif
}

// 打印矢量，一般用于 debug
__device__ __forceinline__ void vct6_print(FLOAT *v) {
#ifdef FLOAT32
    printf("%.15f, %.15f, %.15f, %.15f, %.15f, %.15f\n", v[X], v[Y], v[Z], v[X + DIM], v[Y + DIM], v[Z + DIM]);
#else
    printf("%.15lf, %.15lf, %.15lf, %.15lf, %.15lf, %.15lf\n", v[X], v[Y], v[Z] ,v[X+DIM], v[Y+DIM], v[Z+DIM]);
#endif
}

// 矢量减法
__device__ __forceinline__ void vct_sub(FLOAT *a, FLOAT *b, FLOAT *ret) {
    ret[X] = a[X] - b[X];
    ret[Y] = a[Y] - b[Y];
    ret[Z] = a[Z] - b[Z];
}


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
__device__ __forceinline__ void dB(FLOAT *kl, FLOAT *p0, FLOAT *p, FLOAT *ret) {
    FLOAT r[DIM];
    FLOAT rr;

    vct_sub(p, p0, r); // r = p - p0

    rr = vct_len(r); // rr = abs(r)

    rr = rr * rr * rr; // rr = rr^3

    vct_cross(kl, r, ret); // ret = kl × r

    vct_dot_a_v(1.0 / rr, ret); // ret = (kl × r)/(rr^3)
}

// 计算所有的电流元在 p 点产生的磁场
// number 表示电流元数目
// kls 每 DIM = 3 组表示一个 kl
// p0s 每 DIM = 3 组表示一个 p0
// shared_ret 应该是一个 shared 量，保存返回值
// 调用该方法后，应该同步处理  __syncthreads();
__device__ void current_element_B(FLOAT *kls, FLOAT *p0s, int number, FLOAT *p, FLOAT *shared_ret) {
    int tid = threadIdx.x; // 0-1023 (decide by BLOCK_DIM_X)
    FLOAT db[DIM];
    __shared__ FLOAT s_dbs[DIM * BLOCK_DIM_X];

    vct_zero(s_dbs + tid * DIM);

// 计算每个电流元产生的磁场
    for (int i = tid * DIM; i < number * DIM; i += BLOCK_DIM_X * DIM) {
        dB(kls + i, p0s + i, p, db);
        vct_add_local(s_dbs + tid * DIM, db);
    }

// 规约求和（from https://www.bilibili.com/video/BV15E411x7yT）
    for (int step = BLOCK_DIM_X >> 1; step >= 1; step >>= 1) {
        __syncthreads(); // 求和前同步
        if (tid < step) vct_add_local(s_dbs + tid * DIM, s_dbs + (tid + step) * DIM);
    }

    if (tid == 0) vct_copy(s_dbs, shared_ret);
}


// 计算 QS 在 p 点产生的磁场
// origin xi yi zi 分别是 QS 的局部坐标系
// 这个函数只需要单线程计算
__device__ __forceinline__ void magnet_at_qs(FLOAT *origin, FLOAT *xi, FLOAT *yi, FLOAT *zi,
                                             FLOAT length, FLOAT gradient, FLOAT second_gradient, FLOAT aper_r,
                                             FLOAT *p, FLOAT *ret) {
    FLOAT temp1[DIM];
    FLOAT temp2[DIM];

    vct_sub(p, origin, temp1); // temp1 = p - origin
    temp2[X] = vct_dot_v_v(xi, temp1);
    temp2[Y] = vct_dot_v_v(yi, temp1);
    temp2[Z] = vct_dot_v_v(zi, temp1); // 这时 temp2 就是全局坐标 p 点在 QS 局部坐标系中的坐标

    vct_zero(ret);

    if (temp2[Z] < 0 || temp2[Z] > length) {
        return; // 无磁场
    } else {
        if (
                temp2[X] > aper_r ||
                temp2[X] < -aper_r ||
                temp2[Y] > aper_r ||
                temp2[Y] < -aper_r ||
#ifdef FLOAT32
                sqrtf(temp2[X] * temp2[X] + temp2[Y] * temp2[Y]) > aper_r
#else
                sqrt(temp2[X]*temp2[X]+temp2[Y]*temp2[Y]) > aper_r
#endif
                ) {
            return; // 无磁场
        } else {
            temp1[X] = gradient * temp2[Y] + second_gradient * (temp2[X] * temp2[Y]);
            temp1[Y] = gradient * temp2[X] + 0.5 * second_gradient * (temp2[X] * temp2[X] - temp2[Y] * temp2[Y]);

            vct_dot_a_v_ret(temp1[X], xi, ret);
            vct_dot_a_v_ret(temp1[Y], yi, temp2);
            vct_add_local(ret, temp2);
        }
    }
}


// 整个束线在 p 点产生得磁场（只有一个 QS 磁铁！）
// FLOAT *kls, FLOAT* p0s, int current_element_number 和 CCT 电流元相关
// FLOAT *qs_data 表示 QS 磁铁所有参数，分别是局部坐标系（原点origin,三个轴xi yi zi，长度 梯度 二阶梯度 孔径）
// p 表示要求磁场得全局坐标点
// shared_ret 表示磁场返回值（应该是一个 __shared__）
// 本方法已经完成同步了，不用而外调用 __syncthreads();
__device__ void magnet_with_single_qs(FLOAT *kls, FLOAT *p0s, int current_element_number,
                                      FLOAT *qs_data, FLOAT *p, FLOAT *shared_ret) {
    int tid = threadIdx.x;
    FLOAT qs_magnet[DIM];

    current_element_B(kls, p0s, current_element_number, p, shared_ret);
    __syncthreads(); // 块内同步


    if (tid == 0) {
        // 计算 QS 的磁场确实不能并行
        // 也没有必要让每个线程都重复计算一次
        // 虽然两次同步有点麻烦，但至少只有一个线程束参与运行
        magnet_at_qs(
                qs_data, // origin
                qs_data + 3, //xi
                qs_data + 6, //yi
                qs_data + 9, //zi
                *(qs_data + 12), // len
                *(qs_data + 13), // g
                *(qs_data + 14), // sg
                *(qs_data + 15), // aper r
                p, qs_magnet
        );

        vct_add_local(shared_ret, qs_magnet);
    }
    __syncthreads(); // 块内同步
}

// runge_kutta4 代码和 cctpy 中的 runge_kutta4 一模一样
// Y0 数组长度为 6
// Y0 会发生变化，既是输入也是输出
// 为了分析包络等，会出一个记录全部 YO 的函数
// 这个函数单线程运行

// void (*call)(FLOAT,FLOAT*,FLOAT*) 表示 tn Yn 到 Yn+1 的转义，实际使用中还会带更多参数（C 语言没有闭包）
// 所以这个函数仅仅是原型
__device__ void runge_kutta4(FLOAT t0, FLOAT t_end, FLOAT *Y0, void (*call)(FLOAT, FLOAT *, FLOAT *), FLOAT dt) {
#ifdef FLOAT32
    int number = (int) (ceilf((t_end - t0) / dt));
#else
    int number = (int)(ceil((t_end - t0) / dt));
#endif

    dt = (t_end - t0) / ((FLOAT) (number));
    FLOAT k1[DIM * 2];
    FLOAT k2[DIM * 2];
    FLOAT k3[DIM * 2];
    FLOAT k4[DIM * 2];
    FLOAT temp[DIM * 2];

    for (int ignore = 0; ignore < number; ignore++) {
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


// runge_kutta4_for_magnet_with_single_qs 函数用到的回调
// FLOAT t0, FLOAT* Y0, FLOAT* Y1 微分计算
// 其中 Y = [P, V]
// FLOAT k = particle[E] / particle[RM]; // k: float = particle.e / particle.relativistic_mass
// FLOAT *kls, FLOAT* p0s, int current_element_number, 表示所有电流元
// FLOAT *qs_data 表示一个 QS 磁铁
__device__ void callback_for_runge_kutta4_for_magnet_with_single_qs(
        FLOAT t0, FLOAT *Y0, FLOAT *Y1, FLOAT k,
        FLOAT *kls, FLOAT *p0s, int current_element_number,
        FLOAT *qs_data
) {
    int tid = threadIdx.x;
    __shared__ FLOAT m[DIM]; // 磁场
    magnet_with_single_qs(kls, p0s, current_element_number, qs_data, Y0, m); //Y0 只使用前3项，表示位置。已同步

    if (tid == 0) { // 单线程完成即可
        // ------------ 以下两步计算加速度，写入 Y1 + 3 中 ----------
        // Y0 + 3 是原速度 v
        // Y1 + 3 用于存加速度，即 v × m，还没有乘 k = e/rm
        vct_cross(Y0 + 3, m, Y1 + 3);
        vct_dot_a_v(k, Y1 + 3); // 即 (v × m) * a，并且把积存在 Y1 + 3 中

        // ------------- 以下把原速度复制到 Y1 中 ------------
        vct_copy(Y0 + 3, Y1); // Y0 中后三项，速度。复制到 Y1 的前3项
    }

    __syncthreads(); // 块内同步
}

// 单个粒子跟踪
// runge_kutta4 函数用于 magnet_with_single_qs 的版本，即粒子跟踪
// Y0 即是 [P, v] 粒子位置、粒子速度
// void (*call)(FLOAT,FLOAT*,FLOAT*,FLOAT,FLOAT*,FLOAT*,int,FLOAT*) 改为 callback_for_runge_kutta4_for_magnet_with_single_qs
// 前 3 项 FLOAT,FLOAT*,FLOAT* 和函数原型 runge_kutta4 函数一样，即 t0 Y0 Y1
// 第 4 项，表示 k = particle[E] / particle[RM]; // k: float = particle.e / particle.relativistic_mass
// 第 567 项，FLOAT*,FLOAT*,int 表示所有电流源，FLOAT *kls, FLOAT* p0s, int current_element_number
// 最后一项，表示 qs_data
// particle 表示粒子 (px0, py1, pz2, vx3, vy4, vz5, rm6, e7, speed8, distance9) len = 10
/*__global__*/ __device__ void track_for_magnet_with_single_qs(FLOAT *distance, FLOAT *footstep,
                                                               FLOAT *kls, FLOAT *p0s, int *current_element_number,
                                                               FLOAT *qs_data, FLOAT *particle) {
    int tid = threadIdx.x;
    FLOAT t0 = 0.0; // 开始时间为 0
    FLOAT t_end = (*distance) / particle[SPEED]; // 用时 = 距离/速率

#ifdef FLOAT32
    int number = (int) (ceilf((*distance) / (*footstep)));
#else
    int number = (int)(ceil( (*distance) / (*footstep)));
#endif

    FLOAT dt = (t_end - t0) / ((FLOAT) (number));
    FLOAT k = particle[E] / particle[RM]; // k: float = particle.e / particle.relativistic_mass

    __shared__ FLOAT Y0[DIM * 2]; // Y0 即是 [P, v] 粒子位置、粒子速度，就是 particle 前两项
    __shared__ FLOAT k1[DIM * 2];
    __shared__ FLOAT k2[DIM * 2];
    __shared__ FLOAT k3[DIM * 2];
    __shared__ FLOAT k4[DIM * 2];
    __shared__ FLOAT temp[DIM * 2];

    if (tid == 0) {
        vct6_copy(particle, Y0); // 写 Y0
    }

    for (int ignore = 0; ignore < number; ignore++) {
        __syncthreads(); // 循环前同步

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0, Y0, k1, k, kls, p0s, *current_element_number,
                                                            qs_data); // 已同步


        if (tid == 0) {
            vct6_dot_a_v_ret(dt / 2., k1, temp); // temp = dt / 2 * k1
            vct6_add_local(temp, Y0); // temp =  Y0 + temp
        }
        __syncthreads();

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0 + dt / 2., temp, k2, k, kls, p0s,
                                                            *current_element_number, qs_data);

        if (tid == 0) {
            vct6_dot_a_v_ret(dt / 2., k2, temp); // temp = dt / 2 * k2
            vct6_add_local(temp, Y0); // temp =  Y0 + temp
        }
        __syncthreads();

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0 + dt / 2., temp, k3, k, kls, p0s,
                                                            *current_element_number, qs_data);

        if (tid == 0) {
            vct6_dot_a_v_ret(dt, k3, temp); // temp = dt * k3
            vct6_add_local(temp, Y0); // temp =  Y0 + temp
        }
        __syncthreads();

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0 + dt, temp, k4, k, kls, p0s, *current_element_number,
                                                            qs_data);

        t0 += dt;

        if (tid == 0) {
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

// 写回 particle
    if (tid == 0) {
        vct6_copy(Y0, particle); // 写 Y0
        particle[DISTANCE] = *distance;
    }

    __syncthreads();
}

// 上函数的 global 版本
__global__ void track_for_magnet_with_single_qs_g(FLOAT *distance, FLOAT *footstep,
                                                  FLOAT *kls, FLOAT *p0s, int *current_element_number,
                                                  FLOAT *qs_data, FLOAT *particle) {
    int tid = threadIdx.x;
    FLOAT t0 = 0.0; // 开始时间为 0
    FLOAT t_end = (*distance) / particle[SPEED]; // 用时 = 距离/速率

#ifdef FLOAT32
    int number = (int) (ceilf((*distance) / (*footstep)));
#else
    int number = (int)(ceil( (*distance) / (*footstep)));
#endif

    FLOAT dt = (t_end - t0) / ((FLOAT) (number));
    FLOAT k = particle[E] / particle[RM]; // k: float = particle.e / particle.relativistic_mass

    __shared__ FLOAT Y0[DIM * 2]; // Y0 即是 [P, v] 粒子位置、粒子速度，就是 particle 前两项
    __shared__ FLOAT k1[DIM * 2];
    __shared__ FLOAT k2[DIM * 2];
    __shared__ FLOAT k3[DIM * 2];
    __shared__ FLOAT k4[DIM * 2];
    __shared__ FLOAT temp[DIM * 2];

    if (tid == 0) {
        vct6_copy(particle, Y0); // 写 Y0
    }

    for (int ignore = 0; ignore < number; ignore++) {
        __syncthreads(); // 循环前同步

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0, Y0, k1, k, kls, p0s, *current_element_number,
                                                            qs_data); // 已同步


        if (tid == 0) {
            vct6_dot_a_v_ret(dt / 2., k1, temp); // temp = dt / 2 * k1
            vct6_add_local(temp, Y0); // temp =  Y0 + temp
        }
        __syncthreads();

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0 + dt / 2., temp, k2, k, kls, p0s,
                                                            *current_element_number, qs_data);

        if (tid == 0) {
            vct6_dot_a_v_ret(dt / 2., k2, temp); // temp = dt / 2 * k2
            vct6_add_local(temp, Y0); // temp =  Y0 + temp
        }
        __syncthreads();

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0 + dt / 2., temp, k3, k, kls, p0s,
                                                            *current_element_number, qs_data);

        if (tid == 0) {
            vct6_dot_a_v_ret(dt, k3, temp); // temp = dt * k3
            vct6_add_local(temp, Y0); // temp =  Y0 + temp
        }
        __syncthreads();

        callback_for_runge_kutta4_for_magnet_with_single_qs(t0 + dt, temp, k4, k, kls, p0s, *current_element_number,
                                                            qs_data);

        t0 += dt;

        if (tid == 0) {
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

// 写回 particle
    if (tid == 0) {
        vct6_copy(Y0, particle); // 写 Y0
        particle[DISTANCE] = *distance;
    }

    __syncthreads();
}


__device__ void track_multi_particle_for_magnet_with_single_qs(FLOAT *distance, FLOAT *footstep,
                                                               FLOAT *kls, FLOAT *p0s, int *current_element_number,
                                                               FLOAT *qs_data, FLOAT *particle, int *particle_number) {
    for (int i = 0; i < (*particle_number); i++) {
        track_for_magnet_with_single_qs(distance, footstep, kls, p0s,
                                        current_element_number, qs_data, particle + i * PARTICLE_DIM);
    }
}

__global__ void track_multi_particle_for_magnet_with_single_qs_g(FLOAT *distance, FLOAT *footstep,
                                                                 FLOAT *kls, FLOAT *p0s, int *current_element_number,
                                                                 FLOAT *qs_data, FLOAT *particle,
                                                                 int *particle_number) {
    for (int i = 0; i < (*particle_number); i++) {
        track_for_magnet_with_single_qs(distance, footstep, kls, p0s,
                                        current_element_number, qs_data, particle + i * PARTICLE_DIM);
    }
}


__global__ void track_multi_particle_beamlime_for_magnet_with_single_qs(FLOAT *distance, FLOAT *footstep,
                                                                        FLOAT *kls, FLOAT *p0s,
                                                                        int *current_element_number,
                                                                        FLOAT *qs_data, FLOAT *particle,
                                                                        int *particle_number) {
    int bid = blockIdx.x;
    track_multi_particle_for_magnet_with_single_qs(
            distance, // 全局相同
            footstep, // 全局相同

            kls + MAX_CURRENT_ELEMENT_NUMBER * DIM * bid,
            p0s + MAX_CURRENT_ELEMENT_NUMBER * DIM * bid, // 当前组电流元参数
            current_element_number + bid, // 当前组电流元数目

            qs_data + QS_DATA_LENGTH * bid, // 当前组 QS 参数

            particle + (*particle_number) * PARTICLE_DIM * bid, // 当前组粒子
            particle_number // 全局相同
    );
}