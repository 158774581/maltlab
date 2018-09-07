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
% ABSTRACT:  ���׹켣�滮������ļ�
%
% INPUTS:InterpolationCycle      ���㲽������λs
%        TotalLength             �켣���ȣ�1xN���飬NΪ�켣���ɶ�������λm����rad
%        Vmax                    �켣���ٶ�����ٶȣ�1xN���飬��λm/s ��rad/s
%        Acc                     �����ٶȣ�1xN���飬��λm/s^2 �� rad/s^2
%        Dec                     �����ٶȣ�1xN���飬��λm/s^2 �� rad/s^2
%        AccRam                  ���ٵ�Acc��ʱ�䣬��λs
%        DecRam                  ���ٵ�Dec��ʱ�䣬��λs
%
% OUTPUTS:MotionData             �켣��ʽ����
%         s 	                 �켣�������飬NxM���飬��λm �� rad
%         v 	                 �ٶ����飬NxM���飬��λm/s ��rad/s
%         a 	                 ���ٶ����飬NxM���飬��λm/s^2 �� rad/s^2
%         t 	                 ʱ�����飬Mx1���飬��λs
%         maxtime 	             ����ʱ�䣬��λs

function MotionData = MotionPlanning(TotalLength, Vmax, Accel, Decel, AccelRamp, DecelRamp, SystemSpeed,Ts)

% ȷ����Ҫ������������ؽڿռ䣬���Ӧ�������ѿ����ռ䣬���Ӧλ�á���̬��
AxisNbr = length(TotalLength);

% ��������ʱ�����
for i = 1 : AxisNbr
    TimeSegment(i,:) = GetSection(TotalLength(i), Vmax(i), Accel(i), Decel(i), AccelRamp, DecelRamp,  SystemSpeed,Ts);
end

% Ѱ���˶�ʱ������ἰ��ʱ�����
maxtime = 0;
for i = 1 : AxisNbr
    if TimeSegment(i,8) > maxtime
       MaxTimeSegment =  TimeSegment(i,:);
       MaxtimeAxis = i;
       maxtime = TimeSegment(i,8);
    end
end

% �˴�Ϊ����Эͬ��������˶�ʱ�����¼�����ڵ��˶�����
MotionData.t=MaxTimeSegment;
% MotionData.t=maxtime;
for i = 1 : AxisNbr
    if abs(TotalLength(i)) > 1e-10
%         MotionData.t(i,:) = MaxTimeSegment;
        % ���¼������Ӽ��ٶ�
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

% ABSTRACT:  ����ʱ�����
%   
% INPUTS:   TotalLength   λ��
%           Vmax          ����ٶ�
%           Accel         �����ٶ�
%           Decel         �����ٶ�
%           AccelRamp     �����������ٶ�ʱ��
%           DecelRamp     �����������ٶ�ʱ��
%           InterpolationCycle    ��ɢʱ����
%
% OUTPUTS: 	TimeSegment   �˶�ʱ����� 1x8��������ʾ��ͬ�˶�������ʱ��
%                         
function TimeSegment = GetSection(TotalLength, Vmax, Accel, Decel, AccelRamp, DecelRamp, SystemSpeed,InterpolationCycle)
Ts = InterpolationCycle;
% �ж�λ���Ƿ�Ϊ0
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

% motion profile�����־λ
reverse_flag = false;

% �����ٶε��ܴﵽ������ٶ�С�ڼ��ٶ�ʱ�����з�����
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
    % �������ٶθ���ʧ��λ��
    t11 = (v - a1 * t1) / a1;
    t22 = (v - a2 * t2) / a2;
    S1 = v * t1 + 1/2 * v * t11;
    S2 = v * t2 + 1/2 * v * t22;
    S = S1 + S2;

    % �����ȼ��ٸ���ʧ��λ��
    T22 = (a1 * t1 - a2 * t2) / a2;
    Smid = a1 * t1^2 + a1 * t1 * t2 + 1/2 * (a1 * t1) * T22;

    % �����ȼ��ٸ���ʧ��λ��
    Sshort = a2 * t2^2 + J1 * (sqrt((a2 * t2 * t1) / a1))^3; 

    if d > S    % �ж��Ƿ�������ٶ�
        t3 = (d - S) / v;                  
    else
        if d <= S && d > Smid    % �ж��Ƿ�����ȼ��ٶ�
            [t11,t22] = Cal_Situation_1(d, J1, J2, t1, t2);
            t3 = 0;             
        else
            if d <= Smid && d > Sshort     % �ж��Ƿ�����ȼ��ٶ�   
                [t1, t22] = Cal_Situation_2(d, J1, J2, t1, t2);
                t11 = 0;
                t3 = 0;
            else    % ���١��ȼ��١��ȼ��پ���ʧ
                [t1, t2] = Cal_Situation_3(d, J1, J2);
                t11 = 0;
                t22 = 0;
                t3 = 0;                   
            end
        end
    end
else
    if J1 * (t1^2) >= v && J2 * (t2^2) >= v
        % �������ٶθպ���ʧ��λ��
        t1 = sqrt(v / J1);
        t2 = sqrt(v / J2);
        S = J1 * (t1^3) + J2 * (t2^3);
        if d > S    % �ж��Ƿ�������ٶ�
            t3 = (d - S) / v;
            t11 = 0;
            t22 = 0;
        else     % ���١��ȼ��١��ȼ��پ���ʧ
            [t1, t2] = Cal_Situation_3(d, J1, J2);
            t11 = 0;
            t22 = 0;
            t3 = 0;
        end
    else
        % �������ٶθպ���ʧ��λ��
        t1 = sqrt(v / J1);
        t11 = 0;
        t22 = (v - J2 * t2^2) / a2;
        S1 = v * t1 + 1/2 * v * t11;
        S2 = v * t2 + 1/2 * v * t22;
        S = S1 + S2;
        % �����ȼ��ٶθպ���ʧ��λ��
        T1 = sqrt((a2 * t2 * t1) / a1);
        Sshort = J2 * (t2^3) + J1 * (T1^3);
        if d > S    % �ж��Ƿ�������ٶ�
            t3 = (d - S) / v;
        else
            if d <= S && d >= Sshort   % �ж��Ƿ�����ȼ��ٶ�
                [t1, t22] = Cal_Situation_2(d, J1, J2, t1, t2);          
                t11 = 0;
                t3 = 0;
            else    % ���١��ȼ��١��ȼ��پ���ʧ
                [t1, t2] = Cal_Situation_3(d, J1, J2);
                t11 = 0;
                t22 = 0;
                t3 = 0;
            end
        end       
    end
end

% ������з����������ٴη���
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

% ����ʱ������
k = 1/SystemSpeed;
t1 = t1*k;
t11 = t11*k;
t2 = t2*k;
t22 = t22*k;
t3 = t3*k;

% ÿһ��ʱ��ȡ��
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
    % һԪ���η��������ʽϵ��
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
    % �������˶�ʱ���������ȡ��
%     t(8) = ceil((2*t1+2*t2+t11+t22+t3)/Ts)*Ts;
    % ��ȡ��������ʱ��ƽ�����������ٶ�ʱ��t2
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