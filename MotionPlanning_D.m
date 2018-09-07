% Copyright (c) 2016 by QKM Technology (Dong Guan) Co., Ltd.
%
% +-------------------------------------------------------------------------+
% | The information set forth in this document is the property of			|
% | QKM Technology Co., Ltd and is to be held in trust and confidence		|
% | Publication, duplication, disclosure, or use for any purpose not 		|
% | expressly authorized by QKM in writing is prohibited.					|
% | 																		|
% | The information in this document is subject to change without notice	|
% | and should not be construed as a commitment by QKM.						|
% |																			|
% | QKM makes no warranty as to the suitability of this material			|
% | for use by the recipient, and assumes no responsibility for any			|
% | consequences resulting from such use.									|
% +-------------------------------------------------------------------------+
%		
% Author: Jionhui Lin, 09/2016
%
% ABSTRACT:  三阶轨迹规划最基本文件
%
% INPUTS:InterpolationCycle      计算步长，单位s
%        TotalLength             轨迹长度，1xN数组，N为轨迹自由度数，单位m或者rad
%        Vmax                    轨迹匀速段最大速度，1xN数组，单位m/s 或rad/s
%        Acc                     最大加速度，1xN数组，单位m/s^2 或 rad/s^2
%        Dec                     最大减速度，1xN数组，单位m/s^2 或 rad/s^2
%        AccRam                  加速到Acc的时间，单位s
%        DecRam                  减速到Dec的时间，单位s
%
% OUTPUTS:MotionData             轨迹公式参数
%         s 	                 轨迹长度数组，NxM数组，单位m 或 rad
%         v 	                 速度数组，NxM数组，单位m/s 或rad/s
%         a 	                 加速度数组，NxM数组，单位m/s^2 或 rad/s^2
%         t 	                 时间数组，Mx1数组，单位s
%         maxtime 	             所需时间，单位s

function MotionData = MotionPlanning(TotalLength, Vmax, Accel, Decel, AccelRamp, DecelRamp, SystemSpeed,Ts)

% 确定需要计算轴数（如关节空间，则对应轴数；笛卡尔空间，则对应位置、姿态）
AxisNbr = length(TotalLength);

% 计算各轴的时间参数
for i = 1 : AxisNbr
    TimeSegment(i,:) = GetSection(TotalLength(i), Vmax(i), Accel(i), Decel(i), AccelRamp, DecelRamp,  SystemSpeed,Ts);
end

% 寻找运动时间最长的轴及其时间参数
maxtime = 0;
for i = 1 : AxisNbr
    if TimeSegment(i,8) > maxtime
       MaxTimeSegment =  TimeSegment(i,:);
       MaxtimeAxis = i;
       maxtime = TimeSegment(i,8);
    end
end

% 此处为多轴协同，根据最长运动时间重新计算各节点运动参数
MotionData.t=MaxTimeSegment;
% MotionData.t=maxtime;
for i = 1 : AxisNbr
    if abs(TotalLength(i)) > 1e-10
%         MotionData.t(i,:) = MaxTimeSegment;
        % 重新计算最大加加速度
        [J1, J2] = Cal_ReSolve(TotalLength(i), MotionData.t,Ts);       
        MotionData.j(i,:) = [0 J1 0 -J1 0 -J2 0 J2];
        MotionData.a(i,1) = 0;
        MotionData.v(i,1) = 0;
        MotionData.s(i,1) = 0;
        for k = 1:7
            MotionData.a(i,k+1) = Cal_Accelelation(MotionData.t(k+1) - MotionData.t(k), MotionData.j(i,k+1), MotionData.a(i,k),Ts);
            MotionData.v(i,k+1) = Cal_Velocity(MotionData.t(k+1) - MotionData.t(k), MotionData.j(i,k+1), MotionData.a(i,k), MotionData.v(i,k),Ts);
            MotionData.s(i,k+1) = Cal_Displacement(MotionData.t(k+1) - MotionData.t(k), MotionData.j(i,k+1), MotionData.a(i,k), MotionData.v(i,k), MotionData.s(i,k),Ts);
        end
    else

    MotionData.j(i,:)=zeros(1,8);  
    MotionData.a(i,:)=zeros(1,8);  
    MotionData.v(i,:)=zeros(1,8);  
    MotionData.s(i,:)=zeros(1,8); 
    end
end
end

% ABSTRACT:  计算时间参数
%   
% INPUTS:   TotalLength   位移
%           Vmax          最大速度
%           Accel         最大加速度
%           Decel         最大减速度
%           AccelRamp     上升到最大加速度时间
%           DecelRamp     上升到最大减速度时间
%           InterpolationCycle    离散时间间隔
%
% OUTPUTS: 	TimeSegment   运动时间参数 1x8向量，表示不同运动段所需时间
%                         
function TimeSegment = GetSection(TotalLength, Vmax, Accel, Decel, AccelRamp, DecelRamp, SystemSpeed,InterpolationCycle)
Ts = InterpolationCycle;
% 判断位移是否为0
if abs(TotalLength) < 1e-10
    TimeSegment = zeros(1,8);
    return
end

d = abs(TotalLength);
v = Vmax;
a1 = Accel;
a2 = Decel;
t1 = AccelRamp;
t2 = DecelRamp;

J1 = a1 / t1;
J2 = a2 / t2;

% motion profile反向标志位
reverse_flag = false;

% 当加速段的能达到的最大速度小于减速段时，进行反向处理
if(J1 * t1^2) < (J2 * t2^2)
    reverse_flag = true;
    temp = J1;
    J1 = J2;
    J2 = temp;
    temp = a1;
    a1 = a2;
    a2 = temp;
    temp = t1;
    t1 = t2;
    t2 = temp;
end


if (J1 * t1^2 < v) && (J2 * t2^2 < v)
    % 计算匀速段刚消失的位移
    t11 = (v - a1 * t1) / a1;
    t22 = (v - a2 * t2) / a2;
    S1 = v * t1 + 1/2 * v * t11;
    S2 = v * t2 + 1/2 * v * t22;
    S = S1 + S2;

    % 计算匀加速刚消失的位移
    T22 = (a1 * t1 - a2 * t2) / a2;
    Smid = a1 * t1^2 + a1 * t1 * t2 + 1/2 * (a1 * t1) * T22;

    % 计算匀减速刚消失的位移
    Sshort = a2 * t2^2 + J1 * (sqrt((a2 * t2 * t1) / a1))^3; 

    if d > S    % 判断是否存在匀速段
        t3 = (d - S) / v;                  
    else
        if d <= S && d > Smid    % 判断是否存在匀加速段
            [t11,t22] = Cal_Situation_1(d, J1, J2, t1, t2);
            t3 = 0;             
        else
            if d <= Smid && d > Sshort     % 判断是否存在匀减速段   
                [t1, t22] = Cal_Situation_2(d, J1, J2, t1, t2);
                t11 = 0;
                t3 = 0;
            else    % 匀速、匀加速、匀减速均消失
                [t1, t2] = Cal_Situation_3(d, J1, J2);
                t11 = 0;
                t22 = 0;
                t3 = 0;                   
            end
        end
    end
else
    if J1 * (t1^2) >= v && J2 * (t2^2) >= v
        % 计算匀速段刚好消失的位移
        t1 = sqrt(v / J1);
        t2 = sqrt(v / J2);
        S = J1 * (t1^3) + J2 * (t2^3);
        if d > S    % 判断是否存在匀速段
            t3 = (d - S) / v;
            t11 = 0;
            t22 = 0;
        else     % 匀速、匀加速、匀减速均消失
            [t1, t2] = Cal_Situation_3(d, J1, J2);
            t11 = 0;
            t22 = 0;
            t3 = 0;
        end
    else
        % 计算匀速段刚好消失的位移
        t1 = sqrt(v / J1);
        t11 = 0;
        t22 = (v - J2 * t2^2) / a2;
        S1 = v * t1 + 1/2 * v * t11;
        S2 = v * t2 + 1/2 * v * t22;
        S = S1 + S2;
        % 计算匀减速段刚好消失的位移
        T1 = sqrt((a2 * t2 * t1) / a1);
        Sshort = J2 * (t2^3) + J1 * (T1^3);
        if d > S    % 判断是否存在匀速段
            t3 = (d - S) / v;
        else
            if d <= S && d >= Sshort   % 判断是否存在匀减速段
                [t1, t22] = Cal_Situation_2(d, J1, J2, t1, t2);          
                t11 = 0;
                t3 = 0;
            else    % 匀速、匀加速、匀减速均消失
                [t1, t2] = Cal_Situation_3(d, J1, J2);
                t11 = 0;
                t22 = 0;
                t3 = 0;
            end
        end       
    end
end

% 假如进行反向处理，则需再次反向
if reverse_flag == true
    temp = J1;
    J1 = J2;
    J2 = temp;
    temp = t1;
    t1 = t2;
    t2 = temp;
    temp = t11;
    t11 = t22;
    t22 = temp;
end 

% 进行时间缩放
k = 1/SystemSpeed;
t1 = t1*k;
t11 = t11*k;
t2 = t2*k;
t22 = t22*k;
t3 = t3*k;

% 每一段时间取整
t1 =ceil(t1/Ts);
t11 = ceil(t11/Ts);
t2 = ceil(t2/Ts);
t22 = ceil(t22/Ts);
t3 = ceil(t3/Ts);


TimeSegment = Cal_TimeSection(t1,t11,t2,t22,t3);
  
end

function [t1, t2] = Cal_Situation_3(d, J1, J2)
    t1 = (d/(J1 + J2*(sqrt(J1/J2))^3))^(1/3);
    t2 = sqrt(J1/J2)*t1;
end

function [t1, t22] = Cal_Situation_2(d, J1, J2, t1, t2)
    a2 = J2 * t2;
    A = J1^2;
    B = 2 * J1 * a2;
    C = 2 * J1 * a2 * t2 - J1 * J2 * t2^2;
    D = 0;
    E = - 2 * a2 * d;
    x0 = 0;
    x1 = t1;
    if func(A, B, C, D, E, x0) * func(A, B, C, D, E, x1) > 0
        error('no real root')
    end
    xc = (x0 + x1)/2;
    while abs(func(A, B, C, D, E, xc)) > 0.00001
        if func(A, B, C, D, E, x0) * func(A, B, C, D, E, xc) < 0
            x0 = x0;
            x1 = xc;
        else
            x0 = xc;
            x1 = x1;
        end
        xc = (x0 + x1)/2;
    end
    t1 = xc;
    t22 = (J1*t1^2-J2*t2^2)/a2;
end

function y = func(A, B, C, D, E, xc)
    y = A * xc^4 + B * xc^3 + C * xc^2 + D * xc + E;
end


function [t11,t22] = Cal_Situation_1(d, J1, J2, t1, t2)
    a1 = J1 * t1;
    a2 = J2 * t2;
    % 一元二次方程求根公式系数
    A = a1*(a1 + a2); 
    B = J1 * t1^2 * (2 * a1 + a2) - J2 * t2^2 * a1 + 2 * a1 * a2 * (t1 +t2);
    C = J1 * t1^2 * (2 * a2 * t1 + J1 * t1^2 - J2 * t2^2 + 2 * a2 * t2) - 2 * a2 * d;
    
    deta_b = B ^2 - 4 * A * C; 
    
    if deta_b < 0
        error('no real root')
    end

    t11 = (- B + sqrt(deta_b))/(2*A);
    t22 = (J1*t1^2 + a1*t11 - J2*t2^2)/a2;
end

function t = Cal_TimeSection(t1,t11,t2,t22,t3)
    % 对最终运动时间进行向上取整
%     t(8) = ceil((2*t1+2*t2+t11+t22+t3)/Ts)*Ts;
    % 将取整后多余的时间平均分配给变减速段时间t2
%     t2 = (t(8) - (2*t1+t11+t22+t3))/2;
    t(1) = 0;
    t(2) = t(1) + t1;
    t(3) = t(2) + t11;
    t(4) = t(3) + t1;
    t(5) = t(4) + t3;
    t(6) = t(5) + t2;
    t(7) = t(6) + t22;
    t(8) = t(7) + t2;
end

function a = Cal_Accelelation(t, j, a0,Ts)
    a = a0 + j * t*Ts;
end

function v = Cal_Velocity(t, j, a0, v0,Ts)
    v = v0 + a0 * t*Ts+ 1/2 * j * t^2*Ts*Ts;
end

function s = Cal_Displacement(t, j, a0, v0, s0,Ts)
    s = s0 + v0 * t*Ts+ 1/2 * a0 * t^2*Ts^2 + 1/6 * j * t^3 *Ts^3;
end

function [J1, J2] = Cal_ReSolve(d, t,Ts)
    t1 = t(2) - t(1);
    t11 = t(3) - t(2);
    t3 = t(5) - t(4);
    t2 = t(6) - t(5);
    t22 = t(7) - t(6);
    Vmax = d / (t1 + 1/2 * t11 + t2 + 1/2 * t22 + t3)/Ts;
    J1 = Vmax / (t1^2 + t1 * t11)/Ts^2;
    J2 = Vmax / (t2^2 + t2 * t22)/Ts^2;
end