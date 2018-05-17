%****************************************************************************
%计算3D空间中ECR放电的过程，方法为PIC-MCC
%****************************************************************************

clear variables;
clc;

global m_e q_e m_H q_H
%****************************************************************************
%定义基础参数
%****************************************************************************
m_e = 9.109e-31;        %电子质量
q_e = -1.602e-19;       %电子电量
AMU = 1.661e-27;        %质子质量，单位：kg
EPS0 = 8.854e-12;       %真空中的介电常数，单位：F/m
MU0 = 4*pi()*1e-7;      %真空中的磁导率，单位：H/m
k = 1.38e-23;           %玻尔兹曼常数，单位：J/K
c = 3e8;                %光速，单位：m/s
E_i_Li1 = 5.4;          %Li0->Li1的电离能，单位：eV
E_i_Li2 = 75.77;        %Li1->Li2的电离能，单位：eV
E_i_Li3 = 122.664;      %Li2->Li3的电离能，单位：eV
m_Li = 7*AMU;           %锂单质质量
Z = 0.135;              %Z轴范围，单位：m
R = 0.025;              %放电室半径，单位：m
P = 1;                  %离子源内的压强设置为1Pa
V = pi()*R^2*Z;         %离子源放电室的体积，单位：m^3
T = 273.15+300;         %源内的温度设为300℃
n = P/k/T;              %源内中性原子的数密度，根据理想气体方程P=nkT得出
dt_e = 1.5e-12;         %电子运动时间步长
dt_i = 3.0e-10;         %离子运动时间步长
dg = 0.0005;            %单元网格的长度，单位：m
step_num_i = 600;       %离子的时间步数
N_grid_x = (R-(-R))/0.0005+1;   %x方向的格点数
N_grid_y = (R-(-R))/0.0005+1;   %y方向的格点数
N_grid_z = (Z-0)/0.0005+1;      %z方向的格点数
spwt_e = 1e6;                   %每个宏粒子代表的实际电子个数
spwt_i = 1e6;                   %每个宏粒子代表的实际离子个数
N_e = 1000;                     %初始电子数目为1000个
N_i = 1000;                     %初始离子数目为1个，满足初始电中性条件
max_part = 100000;              %粒子容器的最大值
vth_e = c^2-(c/(-q_e*0.01/m_e/c^2+1)^2);    %电子的热速度，设为2 eV
B_x = zeros(N_grid_x,N_grid_y,N_grid_z);    %格点上的Bx预分配空间
B_y = zeros(N_grid_x,N_grid_y,N_grid_z);    %格点上的By预分配空间
B_z = zeros(N_grid_x,N_grid_y,N_grid_z);    %格点上的Bz预分配空间
E_x = zeros(N_grid_x,N_grid_y,N_grid_z);    %格点上的Ex预分配空间
E_y = zeros(N_grid_x,N_grid_y,N_grid_z);    %格点上的Ey预分配空间
E_z = zeros(N_grid_x,N_grid_y,N_grid_z);    %格点上的Ez预分配空间
B_e = zeros(N_e,3); %电子所在位置处的磁感应强度
B_i = zeros(N_i,3); %离子所在位置处的磁感应强度
E_e = zeros(N_e,3); %电子所在位置处的电场强度
E_i = zeros(N_i,3); %离子所在位置处的电场强度
pos_e = zeros(max_part,3); %给电子的位置预分配空间
vel_e = zeros(max_part,3); %给电子的速度预分配空间
pos_i = zeros(max_part,3); %给离子的位置预分配空间
vel_i = zeros(max_part,3); %给离子的速度预分配空间

load magnetic.txt;         %载入磁场数据
load sigma.txt;            %载入电离截面数据
len_B = length(magnetic);  %确定磁场数据的行数
len_s = length(sigma);     %确定电离截面数据的行数

%将磁场数据都分配到格点上
for i=1:len_B
    in_x=(magnetic(i,1)-(-0.02))/0.0005+1; %x方向网格的格点序号
    in_x=round(in_x); %取整
    in_y=(magnetic(i,2)-(-0.02))/0.0005+1; %y方向网格的格点序号
    in_y=round(in_y); %取整
    in_z=(magnetic(i,3)-0)/0.0005+1; %z方向网格的格点序号
    in_z=round(in_z); %取整
    B_x(in_x,in_y,in_z)=magnetic(i,4); %将x方向的磁感应强度分配到网格上
    B_y(in_x,in_y,in_z)=magnetic(i,5); %将y方向的磁感应强度分配到网格上
    B_z(in_x,in_y,in_z)=magnetic(i,6); %将z方向的磁感应强度分配到网格上
    fprintf('there is %d data left to load, please be patient\n', len_B-i);
end
clc;
fprintf('Congratulations! All the magnetic data has been loaded!\n');


%给定电子的初速度和初位置
vel_e=sampleIsotropicVel(vth_e,N_e); %调用子函数给每个电子分配服从麦克斯韦分布且在4π角度各向同性的速度
theta_pos_e=2*pi()*rand(N_e,1); %选定一个随意角度
R_e=R*rand(N_e,1); %电子的径向位置
pos_e=[R_e.*cos(theta_pos_e),R.*sin(theta_pos_e),Z*rand(N_e,1)+0.011]; %放电室的半径为R，长度为0.042m

%开始主循环
for ts_i=1:step_num_i
    for ts_e=1:dt_i/dt_e %离子运动一步，电子运动dt_i/dt_e步
        
        %首先根据MCC判断是否发生碰撞
        sigma_e=zeros(N_e,1); %电子所在位置处发生电离的截面
        Ek_e=m_e*c^2*(1./sqrt(1-sum(((vel_e/c).^2)'))-1)/(-q_e); %根据相对论动能方程求出所有电子的动能，单位为eV
        Ek_e=Ek_e'; %将动能矩阵进行行列转换，变为N_e行×1列
        
        
        %p=1; %用于遍历每个电子
        
        for p=1:N_e %遍历每个电子
            for cross=1:len_s
                if Ek_e(p,1)<sigma(1,1)
                    sigma_e(p,1)=0;
                    break;
                end
                if sigma(cross,1)<=Ek_e(p,1) && cross<len_s
                    sigma_e(p,1)=1e-6*((sigma(cross+1,2)-sigma(cross,2))/(sigma(cross+1,1)-sigma(cross,1))*(Ek_e(p,1)-sigma(cross,1))+sigma(cross,2)); %通过插值的方式找到电子能量所对应的电离截面
                    P_e=1-exp(-n*sigma_e(p,1)*dt_e); %计算该电子发生电离的概率
                    R_e=rand(); %抽取一个随机数用以和电离概率做比较
                    if R_e<P_e %电离碰撞发生
                        f_s=rand(); %给定一个随机数作为入射电子损耗电离能后分配给散射电子的比例
                        Ek_e_s=(Ek_e(p,1)-E_H)*fs; %碰撞后散射电子的能量
                        Ek_e_n=(Ek_e(p,1)-E_H)*(1-fs); %碰撞后产生的新电子的能量
                        N_e=N_e+1; %发生电离碰撞后，增加一个电子
                        N_i=N_i+1; %发生电离碰撞后，增加一个离子
                        Ek_e(p,1)=Ek_e_s; %将散射后电子的能量赋给原入射电子
                        Ek_e=[Ek_e;Ek_e_n]; %将新产生的电子的能量赋值给最后一个电子
                        v_s2=c^2-(c/(-q_e*Ek_e(p,1)/m_e/c^2+1)^2); %散射电子速度的平方
                        v_n2=c^2-(c/(-q_e*Ek_e(p,1)/m_e/c^2+1)^2); %新产生的电子速度的平方
                        [v_s_2,~]=randfixedsum(3,1,v_s2,0,v_s2); %将散射电子速度的平方分解成x、y、z方向速度的平方和
                        [v_n_2,ignore]=randfixedsum(3,1,v_n2,0,v_n2); %将新产生的电子速度的平方分解为x、y、z方向速度的平方和
                        v_s=randomV(v_s_2); %散射电子在x、y、z方向的分量
                        v_n=randomV(v_n_2); %新产生的电子在x、y、z方向的分量
                        vel_e(p,:)=v_s'; %将散射电子的速度分配给电子速度的矩阵
                        vel_e(N_e,:)=v_n'; %将产生的新电子的速度分配给电子速度的矩阵
                        pos_e(N_e,:)=pos_e(p,:); %新产生的电子的位置与入射电子相同
                        v_i_n=randraw('maxwell',k*T/m_H,1); %生成服从麦克斯韦速率分布的新离子的速率
                        [v_i_v,ignore]=randfixedsum(3,1,v_i_n^2,0,v_i_n^2); %将新生成离子速度的平方分解为x、y、z方向速度的平方
                        v_i=randomV(v_i_v); %将新生成的离子求平方根
                        vel_i(N_i,:)=v_i'; %将新生成的离子的速度分配给电子速度的矩阵
                        pos_i(N_i,:)=pos_e(p,:);%新产生的离子的位置与入射电子相同
                        %p=p+1; %分析下一个电子
                    end
                    break;
                end
                if sigma(cross,1)<=Ek_e(p,1)&&cross==len_s
                    sigma_e(p,1)=1e-6*((0-sigma(cross,2))/(0-sigma(cross,1))*(Ek_e(p,1)-sigma(cross,1))+sigma(cross,2)); %插值找截面
                    P_e=1-exp(-n*sigma_e(p,1)*dt_e); %计算该电子发生电离的概率
                    R_e=rand(); %抽取一个随机数用以和电离概率做比较
                    if R_e<P_e %电离碰撞发生
                        f_s=rand(); %给定一个随机数作为入射电子损耗电离能后分配给散射电子的比例
                        Ek_e_s=(Ek_e(p,1)-E_H)*fs; %碰撞后散射电子的能量
                        Ek_e_n=(Ek_e(p,1)-E_H)*(1-fs); %碰撞后产生的新电子的能量
                        N_e=N_e+1; %发生电离碰撞后，增加一个电子
                        N_i=N_i+1; %发生电离碰撞后，增加一个离子
                        Ek_e(p,1)=Ek_e_s; %将散射后电子的能量赋给原入射电子
                        Ek_e=[Ek_e;Ek_e_n]; %将新产生的电子的能量赋值给最后一个电子
                        v_s2=c^2-(c/(-q_e*Ek_e(p,1)/m_e/c^2+1)^2); %散射电子速度的平方
                        v_n2=c^2-(c/(-q_e*Ek_e(p,1)/m_e/c^2+1)^2); %新产生的电子速度的平方
                        [v_s_2,ignore]=randfixedsum(3,1,v_s2,0,v_s2); %将散射电子速度的平方分解成x、y、z方向速度的平方和
                        [v_n_2,ignore]=randfixedsum(3,1,v_n2,0,v_n2); %将新产生的电子速度的平方分解为x、y、z方向速度的平方和
                        v_s=randomV(v_s_2); %散射电子在x、y、z方向的分量
                        v_n=randomV(v_n_2); %新产生的电子在x、y、z方向的分量
                        vel_e(p,:)=v_s'; %将散射电子的速度分配给电子速度的矩阵
                        vel_e(N_e,:)=v_n'; %将产生的新电子的速度分配给电子速度的矩阵
                        pos_e(N_e,:)=pos_e(p,:); %新产生的电子的位置与入射电子相同
                        v_i_n=randraw('maxwell',k*T/m_H,1); %生成服从麦克斯韦速率分布的新离子的速率
                        [v_i_v,ignore]=randfixedsum(3,1,v_i_n^2,0,v_i_n^2); %将新生成离子速度的平方分解为x、y、z方向速度的平方
                        v_i=randomV(v_i_v); %将新生成的离子求平方根
                        vel_i(N_i,:)=v_i'; %将新生成的离子的速度分配给电子速度的矩阵
                        pos_i(N_i,:)=pos_e(p,:);%新产生的离子的位置与入射电子相同
                    end
                end
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
    end
end
        
        
        
        
        
        
        
        
