% code created by       陈鹤鸣
% code re-factoring by  秦斌
% version 0.9, 2020-April-8

% TODO: 所有标量改成矢量！

close all;clear;clc


%-----------输入弯转CCT的设计参数-------------
global n; 
global N;
global I1;
global I2;
global ln;
global wn;
global step;

global miu;

n=1;                    %输入所设计的弯转CCT极数，n=1对应二极，n=2对应四极，以此类推
Rho=1.5;                  %CCT Bending Curvature, unit: m
phi0=0.01;              %CCT的匝间距 phi0, unit: m

disp('Axial Turn Separation (mm): ')
disp(phi0*Rho*1000)

r1=0.06;                %内层CCT/Layer1 孔径半径, unit: m
r2=0.08;                %外层CCT/Layer2 孔径半径, unit: m
alpha1=20.0;            %内层CCT/Layer1 倾斜角 Alpha1, unit: deg.
alpha2=180-alpha1;      %外层CCT/Layer1 倾斜角 Alpha2=180-Alpha1, unit: deg.
alpha1=alpha1*pi/180;   %将倾斜角由角度制转化为弧度制
alpha2=alpha2*pi/180;

N = 5;                  %CCT的匝数
I1 = 12000;             %CCT/Layer1, Current I1
I2 = -I1;               %CCT/Layer2, Current I2; Default I2=-I1

%-----------多线模型的设计参数-------------
length=0.006;            %线圈横截面的长度(m)
width=0.004;             %线圈横截面的宽度(m)
ln=1;                   %线圈横截面长度方向线的数量
wn=1;                   %线圈横截面宽度方向线的数量

delta_w=width/2/wn;     %线圈横截面宽度方向上两线间距离的一半
delta_l=length/2/ln;    %线圈横截面长度方向上两线间距离的一半
I1=I1/ln/wn;            %多线模型每根线上分摊的电流值
I2=I2/ln/wn;

%-----------圆环坐标系下弯转CCT路径方程参数-------------
a1=sqrt(Rho*Rho-r1*r1);                 % 内层弯转CCT双极坐标系的极点
eta01=asinh(sqrt((Rho/r1)^2-1));      % 内层弯转CCT双极坐标系η0
Cn1_0=cot(alpha1)/n/sinh(eta01);    % 内层弯转CCT对应的2n极磁场系数

a2=sqrt(Rho*Rho-r2*r2);                 % 外层弯转CCT双极坐标系的极点
eta02=asinh(sqrt((Rho/r2)^2-1));      % 外层弯转CCT双极坐标系η0
Cn2_0=cot(alpha2)/n/sinh(eta02);    % 外层弯转CCT对应的2n极磁场系数
step=pi/180;                        %  CCT的步距
miu=4*pi*10^(-7);

total_steps = int32(2.0*pi*N/step);
disp('Total Path Points: ');
disp(total_steps);


%% Two Layers Curved CCT Path with Multiple lines Generation!
PathCenter = zeros(total_steps+1,3);            % define array to contain all center points in CCT path
PathCenterVector_t = zeros(total_steps+1,3);      % define array to contain all vector (unit) in CCT central path; -->t
PathCenterVector_r = zeros(total_steps+1,3);      % define array to contain all vector (unit) in CCT central path; -->r
PathCenterVector_b = zeros(total_steps+1,3);      % define array to contain all vector (unit) in CCT central path; -->b

MultiLines_1 = zeros(total_steps+1,ln,wn,3);      % define array to contain multiple lines (ln*wn) in CCT path
MultiLines_2 = zeros(total_steps+1,ln,wn,3);      % define array to contain multiple lines (ln*wn) in CCT path

PathCenter(1,:) = [a1*sinh(eta01)/(cosh(eta01)-cos(0))*cos(Cn1_0*sin(1*0)+phi0/2/pi*0)
              a1*sinh(eta01)/(cosh(eta01)-cos(0))*sin(Cn1_0*sin(1*0)+phi0/2/pi*0)
              a1*sin(0)/(cosh(eta01)-cos(0))];

for ksi=step:step:N*2*pi %弯转CCT是以变量ξ细分的
    phi1 = Cn1_0*sin(n*ksi)+phi0/2/pi*ksi;  % 弯转CCT的变量 φ(ξ)
    k1 = cosh(eta01)-cos(ksi);  %  弯转CCT路径方程中的辅助量 k(ξ)
    i=int32(ksi/step);
    
    PathCenter(i+1,:) = [a1*sinh(eta01)/k1*cos(phi1)
        a1*sinh(eta01)/k1*sin(phi1)
        a1*sin(ksi)/k1];
    
    PathCenterVector_t(i,:) = PathCenter(i+1,:) - PathCenter(i,:);
    PathCenterVector_t(i,:) = PathCenterVector_t(i,:) / norm(PathCenterVector_t(i,:));

    [angle, rho] = cart2pol(PathCenter(i,1), PathCenter(i,2));

    cur_CCT_Center = [Rho*cos(angle) Rho*sin(angle) 0];
    PathCenterVector_r(i,:) = PathCenter(i,:) - cur_CCT_Center;
    PathCenterVector_r(i,:) = PathCenterVector_r(i,:) / norm(PathCenterVector_r(i,:));
    
    PathCenterVector_b(i,:)=cross(PathCenterVector_t(i,:), PathCenterVector_r(i,:)); %内层弯转CCT第i个点的副法向矢量
    PathCenterVector_b(i,:) = PathCenterVector_b(i,:) / norm(PathCenterVector_b(i,:));
    
    for m=1:ln
        delta_r1=(2*(m-ln/2)-1)*delta_l; %多线模型每根线在径向上的坐标
        for k=1:wn
            delta_p1=(2*(k-wn/2)-1)*delta_w; %多线模型每根线在副法向上的坐标
            
            MultiLines_1(1,m,k,:) = [a1*sinh(eta01)/(cosh(eta01)-cos(0))*cos(Cn1_0*sin(1*0)+phi0/2/pi*0) + delta_r1*PathCenterVector_r(1,1) + delta_p1*PathCenterVector_b(1,1)
                a1*sinh(eta01)/(cosh(eta01)-cos(0))*sin(Cn1_0*sin(1*0)+phi0/2/pi*0) + delta_r1*PathCenterVector_r(1,2) + delta_p1*PathCenterVector_b(1,2)
                a1*sin(0)/(cosh(eta01)-cos(0)) + delta_r1*PathCenterVector_r(1,3) + delta_p1*PathCenterVector_b(1,3)
                ];
            
            MultiLines_1(i+1,m,k,:) = [a1*sinh(eta01)/k1*cos(phi1) + delta_r1*PathCenterVector_r(i,1) + delta_p1*PathCenterVector_b(i,1)
                a1*sinh(eta01)/k1*sin(phi1) + delta_r1*PathCenterVector_r(i,2) + delta_p1*PathCenterVector_b(i,2)
                a1*sin(ksi)/k1 + delta_r1*PathCenterVector_r(i,3) + delta_p1*PathCenterVector_b(i,3)
                ];
        
         
        end
    end
    
    
end

figure(1)
for m=1:ln
    for k=1:wn
        plot3(MultiLines_1(:,m,k,1), MultiLines_1(:,m,k,2), MultiLines_1(:,m,k,3));
        hold on
    end
end

axis equal




for ksi=step:step:N*2*pi %弯转CCT是以变量ξ细分的
    phi2 = Cn2_0*sin(n*ksi)+phi0/2/pi*ksi;  % 弯转CCT的变量 φ(ξ)
    k2 = cosh(eta02)-cos(ksi);  %  弯转CCT路径方程中的辅助量 k(ξ)
    i=int32(ksi/step);
    
    PathCenter(i+1,:) = [a2*sinh(eta02)/k2*cos(phi2)
        a2*sinh(eta02)/k1*sin(phi2)
        a2*sin(ksi)/k2];
    
    PathCenterVector_t(i,:) = PathCenter(i+1,:) - PathCenter(i,:);
    PathCenterVector_t(i,:) = PathCenterVector_t(i,:) / norm(PathCenterVector_t(i,:));

    [angle, rho] = cart2pol(PathCenter(i,1), PathCenter(i,2));
    
    cur_CCT_Center = [Rho*cos(angle) Rho*sin(angle) 0];
    PathCenterVector_r(i,:) = PathCenter(i,:) - cur_CCT_Center;
    PathCenterVector_r(i,:) = PathCenterVector_r(i,:) / norm(PathCenterVector_r(i,:));

    PathCenterVector_b(i,:)=cross(PathCenterVector_t(i,:), PathCenterVector_r(i,:)); %内层弯转CCT第i个点的副法向矢量
    PathCenterVector_b(i,:) = PathCenterVector_b(i,:) / norm(PathCenterVector_b(i,:));

    for m=1:ln
        delta_r2=(2*(m-ln/2)-1)*delta_l; %多线模型每根线在径向上的坐标
        for k=1:wn
            delta_p2=(2*(k-wn/2)-1)*delta_w; %多线模型每根线在副法向上的坐标
            
            MultiLines_2(1,m,k,:) = [a2*sinh(eta02)/(cosh(eta02)-cos(0))*cos(Cn2_0*sin(1*0)+phi0/2/pi*0) + delta_r2*PathCenterVector_r(1,1) + delta_p2*PathCenterVector_b(1,1)
                a2*sinh(eta02)/(cosh(eta02)-cos(0))*sin(Cn2_0*sin(1*0)+phi0/2/pi*0) + delta_r2*PathCenterVector_r(1,2) + delta_p2*PathCenterVector_b(1,2)
                a2*sin(0)/(cosh(eta02)-cos(0)) + delta_r2*PathCenterVector_r(1,3) + delta_p2*PathCenterVector_b(1,3)
                ];
            
            MultiLines_2(i+1,m,k,:) = [a2*sinh(eta02)/k2*cos(phi2) + delta_r2*PathCenterVector_r(i,1) + delta_p2*PathCenterVector_b(i,1)
                a2*sinh(eta02)/k2*sin(phi2) + delta_r2*PathCenterVector_r(i,2) + delta_p2*PathCenterVector_b(i,2)
                a2*sin(ksi)/k2 + delta_r2*PathCenterVector_r(i,3) + delta_p2*PathCenterVector_b(i,3)
                ];
         
        end
    end
end

figure(1)
for m=1:ln
    for k=1:wn
        plot3(MultiLines_2(:,m,k,1), MultiLines_2(:,m,k,2), MultiLines_2(:,m,k,3));
        hold on
    end
end

axis equal
% PAUSE!!!!

%% B field of Two Layers Curved CCT Calculation!
    Rho_calculation=1.5;              %输入所计算的圆弧的半径
    Z=0;                   %输入所计算的圆弧在笛卡尔坐标系下的Z值
    start_angle=-10.0;       %计算的初始点角度(°)
    end_angle=15.0;         %计算的结束点角度(°)
    point_num=26;           %所需计算圆弧上的节点数目
    angle=linspace(start_angle, end_angle, point_num);
    
    %CalcuPoint
    
    delta_angle=(end_angle-start_angle)/(point_num-1);
    line_step = Rho_calculation * delta_angle * pi/180.0;
    dB1 = zeros(total_steps+1,3);
    sumdB1 = zeros(ln,wn,3);
    B1 = zeros(point_num,3);
    
    dB2 = zeros(total_steps+1,3);
    sumdB2 = zeros(ln,wn,3);
    B2 = zeros(point_num,3);
    
    B_total = zeros(point_num,3);

    for t=1:point_num
        CalcuPoint=[Rho_calculation*cos((start_angle+delta_angle*(t-1))/180*pi); Rho_calculation*sin((start_angle+delta_angle*(t-1))/180*pi); Z];
        %disp(CalcuPoint)
        %disp(MultiLines_1(i,m,k,:))
            
        for m=1:ln
            for k=1:wn
                for ksi=step:step:N*2*pi %弯转CCT是以变量ξ细分的
                    i=int32(ksi/step);
                    
                    %%Layer #1
                    Current_vector = [MultiLines_1(i+1,m,k,1) - MultiLines_1(i,m,k,1)
                        MultiLines_1(i+1,m,k,2) - MultiLines_1(i,m,k,2)
                        MultiLines_1(i+1,m,k,3) - MultiLines_1(i,m,k,3)
                        ];
                    
                    Current_Center_Point = [MultiLines_1(i+1,m,k,1) + MultiLines_1(i,m,k,1)
                        MultiLines_1(i+1,m,k,2) + MultiLines_1(i,m,k,2)
                        MultiLines_1(i+1,m,k,3) + MultiLines_1(i,m,k,3)
                        ]/2.0;
                    
%                     Current_Center_Point = [MultiLines_1(i,m,k,1)
%                         MultiLines_1(i,m,k,2)
%                         MultiLines_1(i,m,k,3)
%                         ];
                    
                    CurrentToPoint_vector = CalcuPoint - Current_Center_Point;
                    radial_vector_norm = norm(CurrentToPoint_vector);
                    er=CurrentToPoint_vector / radial_vector_norm;
                    multi=cross(Current_vector, er);
                    dB1(i,:)=miu/4/pi*I1*multi/radial_vector_norm/radial_vector_norm; %Biot-Savart计算每个微分的磁场

                    
                    %% LAYER #2
                    % Current_vector 采用此种方式，计算的是 (i+1) - (i)个点的矢量，不完全== 第i个点的矢量
                    Current_vector = [MultiLines_2(i+1,m,k,1) - MultiLines_2(i,m,k,1)
                        MultiLines_2(i+1,m,k,2) - MultiLines_2(i,m,k,2)
                        MultiLines_2(i+1,m,k,3) - MultiLines_2(i,m,k,3)
                        ];
                    
                    Current_Center_Point = [MultiLines_2(i+1,m,k,1) + MultiLines_2(i,m,k,1)
                        MultiLines_2(i+1,m,k,2) + MultiLines_2(i,m,k,2)
                        MultiLines_2(i+1,m,k,3) + MultiLines_2(i,m,k,3)
                        ]/2.0;
                    
                    %% TEST, Use Chen's Method for check
%                     Current_Center_Point = [MultiLines_2(i,m,k,1)
%                         MultiLines_2(i,m,k,2)
%                         MultiLines_2(i,m,k,3)
%                         ];
                    

                    CurrentToPoint_vector = CalcuPoint - Current_Center_Point;
                    radial_vector_norm = norm(CurrentToPoint_vector);
                    er=CurrentToPoint_vector / radial_vector_norm;
                    multi=cross(Current_vector, er);
                    dB2(i,:)=miu/4/pi*I2*multi/radial_vector_norm/radial_vector_norm; %Biot-Savart计算每个微分的磁场
                
                end
            sumdB1(k,m,:)=sum(dB1,1);
            sumdB2(k,m,:)=sum(dB2,1);
            end
        end
        B1(t,:) = [sum(sum(sumdB1(:,:,1)))  sum(sum(sumdB1(:,:,2)))   sum(sum(sumdB1(:,:,3)))];
        B2(t,:) = [sum(sum(sumdB2(:,:,1)))  sum(sum(sumdB2(:,:,2)))   sum(sum(sumdB2(:,:,3)))];
        %B_total(t,:) =  B1(t,:) +  B2(t,:);
        B_total(t,:) = [B1(t,1)+B2(t,1)  B1(t,2)+B2(t,2)   B1(t,3)+B2(t,3)];
        %B(t,:) = sum(sum(sumdB1));
    end

    disp(B_total)
    %%Print Field Integrals for Bx, By, Bz!
    disp('Field integrals for Bx, By, Bz [unit: T.m]:')
    disp(sum(B_total) * line_step)
    
    figure(2)
    hold on

    line_Bx = plot(angle, B_total(:,1),'-^g','LineWidth',1);
    line_By = plot(angle, B_total(:,2),'-xr','LineWidth',2);
    line_Bz = plot(angle, B_total(:,3),'--ob', 'LineWidth',2);
    
    xlabel('s (m)','Fontname', 'Times New Roman','FontSize',18);
    ylabel('B (T)','Fontname', 'Times New Roman','FontSize',18);
    legend('Bx', 'By', 'Bz')
    set(gca,'FontSize',13);
    grid on
    

