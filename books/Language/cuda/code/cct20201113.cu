#include <stdio.h>
#include <math.h>

#include <cuda_runtime.h>

#include <helper_cuda.h>

#define MM 0.001f
#define DIM 3
#define PI 3.1415927f
#define X 0
#define Y 1
#define Z 2
#define Proton_Charge_Quantity 1.6021766208e-19f
#define Proton_Static_MassKg 1.672621898e-27f
#define Proton_Static_MassMeV 938.2720813f
#define Light_Speed 299792458.0f
#define RUN_STEP 0.001f

__device__ __forceinline__ void add3d(float* a, float* b, float* ret)
{
    ret[X] = a[X] + b[X];
    ret[Y] = a[Y] + b[Y];
    ret[Z] = a[Z] + b[Z];
}

__device__ __forceinline__ void add3d_local(float* a_local, float* b)
{
    a_local[X] += b[X];
    a_local[Y] += b[Y];
    a_local[Z] += b[Z];
}

__device__ __forceinline__ void sub3d(float* a, float* b, float* ret)
{
    ret[X] = a[X] - b[X];
    ret[Y] = a[Y] - b[Y];
    ret[Z] = a[Z] - b[Z];
}

__device__ __forceinline__ void copy3d(float* src, float* des)
{
    des[X] = src[X];
    des[Y] = src[Y];
    des[Z] = src[Z];
}

__device__ __forceinline__ void cross3d(float* a, float* b, float* ret)
{
    ret[X] = a[Y] * b[Z] - a[Z] * b[Y];
    ret[Y] = -a[X] * b[Z] + a[Z] * b[X];
    ret[Z] = a[X] * b[Y] - a[Y] * b[X];
}

__device__ __forceinline__ void dot_a_v(float a, float* v)
{
    v[X] *= a;
    v[Y] *= a;
    v[Z] *= a;
}

__device__ __forceinline__ void dot_a_v_ret(float a, float* v, float* ret)
{
    ret[X] = v[X] * a;
    ret[Y] = v[Y] * a;
    ret[Z] = v[Z] * a;
}

__device__ __forceinline__ float dot_v_v(float* v1, float* v2)
{
    return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}

__device__ __forceinline__ float len3d(float* v)
{
    return sqrtf(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
}

__device__ __forceinline__ void neg3d(float* v)
{
    v[X] = -v[X];
    v[Y] = -v[Y];
    v[Z] = -v[Z];
}

// 注意，这里计算的不是电流元的磁场，还需要乘以 电流 和 μ0/4π (=1e-7)
// 2020年11月11日 测试通过
__device__  void dB(float* p0, float* p1, float* p, float* ret)
{
    float p01[DIM];
    float r[DIM];
    float rr;

    sub3d(p1, p0, p01); // p01 = p1 - p0

    add3d(p0, p1, r); // r = p0 + p1

    dot_a_v(0.5, r); // r = (p0 + p1)/2

    sub3d(p, r, r); // r = p - r

    rr = len3d(r); // rr = len(r)

    cross3d(p01, r, ret); // ret = p01 x r

    rr = 1.0 / rr / rr / rr; // changed

    dot_a_v(rr, ret); // rr . (p01 x r)
}

// line = float[length][3]。计算导线 line 在 p 点产生的磁场，length 表示 line 中点数目，返回值放到 ret
__device__ void magnet_at_point(float** line, int length, float current, float* p, float* ret)
{
    int i;
    float db[3];

    ret[X] = 0.0f;
    ret[Y] = 0.0f;
    ret[Z] = 0.0f;

    for (i = 0; i < length - 1; i++)
    {
        dB(line[i], line[i + 1], p, db);
        add3d_local(ret, db);
    }

    dot_a_v(current * 1e-7, ret);
}

// 粒子走一步 m 磁场，p 位置，v 速度，rm 动质量，sp 速率
 __device__  __forceinline__  void particle_run_step(float* m, float* p, float* v, float run_mass, float speed)
{
    float a[3]; // 加速度
    float t;    // 运动时间
    float d[3]; // 位置变化 速度变化

    // q v b
    cross3d(v, m, a); // a = v*b

    dot_a_v(Proton_Charge_Quantity / run_mass, a); // a = q v b / mass 加速度

    t = RUN_STEP / speed; // 运动时长

    dot_a_v_ret(t, v, d); // d = t v 位置变化

    add3d_local(p, d); // p+=d

    dot_a_v_ret(t, a, d); // d = t a 速度变化

    add3d_local(v, d); // v+=d
}

// 粒子在单导线磁场中运动，用于测试
// len 运动距离，p 位置，v 速度，rm 动质量，sp 速率，line[des_len][3] 导线，des_len 导线点数，current 电流
 __device__    void particle_run_len_one_line(float len, float* p, float* v, float run_mass, float speed, float** line, int des_len, float current)
{
    float distance = 0.0f;
    float m[3]; // 磁场
    while (distance < len)
    {
        magnet_at_point(line, des_len, current, p, m);
        particle_run_step(m, p, v, run_mass, speed);
        distance += RUN_STEP;
    }
}

 __device__  void linspace(float start, float end, int number, float* ret)
{
    int i;
    float d = (end - start) / (number - 1);

    for (i = 0; i < number - 1; i++)
    {
        ret[i] = start + d * i;
    }

    // 最后一个数赋值，减少计算误差
    ret[number - 1] = end;
}

 // 不可以这样分配显存
 __global__    void test_solenoid() 
{
    int i;
    float t;
    float r = 0.01f;
    float length = 0.1f;
    int n = 20;
    int pas = 360;
    float total_theta = n * 2 * PI;

    float current = 10000.0f;

    int des_len = n * pas;

    float** lines = (float**)malloc(des_len * sizeof(float*));

    float* lins = (float*)malloc(des_len * sizeof(float));

    linspace(0, total_theta, des_len, lins);

    for (i = 0; i < des_len; i++)
    {
        lines[i] = (float*)malloc(DIM * sizeof(float));
        t = lins[i];
        lines[i][X] = cosf(t) * r;
        lines[i][Y] = sinf(t) * r;
        lines[i][Z] = t / total_theta * length;
    }

    float p[3] = { 0, 0, 0 };
    float v[3] = { 0.0, 0.0, 1.839551780274753E8 };

    particle_run_len_one_line(10.0, p, v, 2.1182873748205775E-27, 1.839551780274753E8, lines, des_len, current);

    printf("px=%f --cuda  ", p[X]);
    printf("py=%f --cuda ", p[Y]);
    printf("py=%f --cuda  ", p[Z]);
    printf("vx=%f --cuda  ", v[X]);
    printf("vy=%f --cuda  ", v[Y]);
    printf("vz=%f --cuda  ", v[Z]);

    free(lins);
    for (i = 0; i < des_len; i++)
    {
        free(lines[i]);
    }

    free(lines);
}

 // 测试函数，点 p 在螺线管磁场
 //part_deg 1度 2度 3度，必须是360约数
 __device__ void magnet_at_solenoid(float current,float small_r, float length, float number_wind, int part_deg, float* p,float* sin_table, float* ret) {
     int total_deg = 360 * number_wind;
     float m[3];
     
     int per_deg = 0;
     int cur_deg = part_deg;

     float pre_p[3];
     pre_p[X] = small_r * sin_table[90 + per_deg]; // sin(90+a)=cos(a)
     pre_p[Y] = small_r * sin_table[per_deg];
     pre_p[Z] = length * per_deg / 360.0 / number_wind;


     float cur_p[3];

     ret[X] = 0.0f;
     ret[Y] = 0.0f;
     ret[Z] = 0.0f;

     while (cur_deg <= total_deg) {
         cur_p[X] = small_r * sin_table[(90 + cur_deg)%360];
         cur_p[Y] = small_r * sin_table[cur_deg%360];
         cur_p[Z] = length * cur_deg / 360.0 / number_wind;

         // (float* p0, float* p1, float* p, float* ret)
         dB(pre_p, cur_p, p, m);
         
         add3d_local(ret, m);

         copy3d(cur_p, pre_p);

         cur_deg += part_deg;
     }

     dot_a_v(current * 1e-7, ret);
     
 }

 __global__ void magnet_at_solenoid_test_p000(float* sin_table) {

     float p[3] = { 0,0,0 };
     float m[3];
    
     magnet_at_solenoid(10000.0f, 0.01f, 0.1f, 20, 1, p, sin_table, m);

     printf("mx=%f --cuda  ", m[X]);
     printf("my=%f --cuda  ", m[Y]);
     printf("mz=%f --cuda  ", m[Z]);
 }

 __global__    void test_solenoid_no_malloc(float* sin_table)
 {
     float r = 0.01f;
     float length = 0.1f;
     int n = 20;
     float total_theta = n * 2 * PI;

     float current = 10000.0f;


     float run_len = 1.0f;

     float p[3] = { 0, 0, 0 };
     float v[3] = { 0.0, 0.0, 1.839551780274753E8 };
     float rm = 2.1182873748205775E-27;
     float speed = 1.839551780274753E8;
     float m[3];

     printf("px=%f --cuda  \n", p[X]);
     printf("py=%f --cuda \n", p[Y]);
     printf("py=%f --cuda  \n", p[Z]);
     printf("vx=%f --cuda  \n", v[X]);
     printf("vy=%f --cuda  \n", v[Y]);
     printf("vz=%f --cuda  \n", v[Z]);

     

     float distance = 0.0f;

     while (distance < run_len) {
         magnet_at_solenoid(current, r, length, n, 1, p, sin_table, m);
         particle_run_step(m, p, v, rm, speed);
         distance += RUN_STEP;
     }


     printf("px=%f --cuda  \n", p[X]);
     printf("py=%f --cuda \n", p[Y]);
     printf("py=%f --cuda  \n", p[Z]);
     printf("vx=%f --cuda  \n", v[X]);
     printf("vy=%f --cuda  \n", v[Y]);
     printf("vz=%f --cuda  \n", v[Z]);
 }




/**
 * Host main routine
 */
int
main(void)
{
    int i;
    float sin_table[360];
    for (i = 0; i < 360; i++) {
        sin_table[i] = sin(((float)i) / 180.0f * PI);
    }

    float* d_sin_tb;

    cudaMalloc((void**)&d_sin_tb, 360 * sizeof(float));

    cudaMemcpy(d_sin_tb, sin_table, 360 * sizeof(float), cudaMemcpyHostToDevice);

    test_solenoid_no_malloc <<<256, 24>>>(d_sin_tb);

    cudaDeviceSynchronize();

    cudaFree(d_sin_tb);

    printf("Done\n");
    return 0;
}

