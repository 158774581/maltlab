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
% DESCRIPTION:用牛顿-欧拉迭代动力学算法，由关节运动计算关节力矩。第一部分是对每个
%             连杆应用牛顿-欧拉方程，从连杆1到连杆n向外迭代计算连杆的速度和加速度。
%             第二部分是从连杆n到连杆1向内迭代计算连杆间的相互作用力和力矩以及关节
%             驱动力矩。
%
% AUTHOR : Hui Du,09/2016
%
% ABSTRACT：这是Helios六轴机器人逆动力学函数，通用函数，适用于串联机器人
%           用牛顿-欧拉迭代动力学算法，由关节运动计算关节力矩。第一部分是对每个
%           连杆应用牛顿-欧拉方程，从连杆1到连杆n向外迭代计算连杆的速度和加速度。
%           第二部分是从连杆n到连杆1向内迭代计算连杆间的相互作用力和力矩以及关节
%           驱动力矩。
% 
% INPUT：q                关节角度，1xN矩阵，单位rad
%        qd               关节角速度，1xN，单位rad/s
%        qdd              关节角加速度，1xN矩阵，单位rad/(s^2)
%        config           机器人配置结构体
%
% OUTPUT: Tau             关节力矩矩阵，1xN，单位Nm       
% 
% REFERENCE：《机器人学导论》 
%
function Tau = Idyn( q,qd,qdd,config)
% setting 
Links=config.Links;
n=config.type;
% 1-n个连杆的质量特性
for i=1:n-1
     m(i)=Links(i).m;
    if Links(i).sigma == 1
        Links(i).d = q(i);
        q(i) = 0; %%%需验证这一步
    end
    if i==1
        P(:,i)=[0;0;0];
        R(:,:,i)=rotz(q(i));
    else
        P(:,i)=[Links(i-1).a;0;Links(i-1).d];
        R(:,:,i)=tr2rt(Transformation(-Links(i-1).offset,Links(i-1)))*rotz(q(i)+Links(i).offset);
    end
    R0(:,:,i)=tr2rt(Transformation(-Links(i).offset,Links(i)))*rotz(Links(i+1).offset);
    I(:,:,i)=R0(:,:,i)*Links(i).I/R0(:,:,i);
    Pc(:,i)=R0(:,:,i)*Links(i).r'+rotz(Links(i+1).offset)\[Links(i).a;0;Links(i).d];
end

%末端执行器的质量特性，加负载
m(n)=config.Load.m;
Pc(:,n)=config.Load.r;
P(:,n)=[0;0;0];

I(:,:,n)=config.Load.I;
R(:,:,n)=tr2rt(Transformation(-Links(n-1).offset,Links(n-1)))*rotz(q(n)+Links(n).offset);

% m(n)=Links(n).m;
% Pc(:,n)=Links(n).r;
% P(:,n)=[0;0;0];
% I(:,:,n)=Links(n).I;
% R(:,:,n)=tr2rt(Transformation(Links(n-1).offset,Links(n-1)))*rotz(q(n)+Links(n).offset);
%%%%%%%%%%
R(:,:,n+1)=eye(3);
P(:,n+1)=[0;0;0];

Z=[0;0;1];
% 外推
w0=[0;0;0];
dw0=[0;0;0];
dv0=[0;0;9.81];
% dv0=[0;0;0];
for i=1:n
    if i==1
        if Links(i).sigma == 0
            w(:,i)=R(:,:,i)'*w0+qd(i)*Z;
            dw(:,i)=R(:,:,i)'*dw0+cross(R(:,:,i)'*w0,qd(i)*Z)+qdd(i)*Z;
            dv(:,i)=R(:,:,i)'*(cross(dw0,P(:,i))+cross(w0,cross(w0,P(:,i)))+dv0);
        else
            w(:,i) = R(:,:,i)'*w0;
            dw(:,i) = R(:,:,i)'* dw0;
            dv(:,i) = R(:,:,i)'*(cross(dw0,P(:,i))+cross(w0,cross(w0,P(:,i)))+dv0)+...
                2*cross(w(:,i),qd(i)*Z) + qdd(i)*Z;
        end
    else
        if Links(i).sigma == 0
            w(:,i)=R(:,:,i)'*w(:,i-1)+qd(i)*Z;
            dw(:,i)=R(:,:,i)'*dw(:,i-1)+cross(R(:,:,i)'*w(:,i-1),qd(i)*Z)+qdd(i)*Z;
            dv(:,i)=R(:,:,i)'*(cross(dw(:,i-1),P(:,i))+cross(w(:,i-1),cross(w(:,i-1),P(:,i)))+dv(:,i-1));
        else
            w(:,i) = R(:,:,i)'*w(:,i-1);
            dw(:,i) = R(:,:,i)'* dw(:,i-1);
            dv(:,i) = R(:,:,i)'*(cross(dw(:,i-1),P(:,i))+cross(w(:,i-1),cross(w(:,i-1),P(:,i)))+dv(:,i-1))+...
                2*cross(w(:,i),qd(i)*Z) + qdd(i)*Z;
        end
    end
    dvc(:,i)=cross(dw(:,i),Pc(:,i))+cross(w(:,i),cross(w(:,i),Pc(:,i)))+dv(:,i);
    F(:,i)=m(i)*dvc(:,i);
    N(:,i)=I(:,:,i)*dw(:,i)+cross(w(:,i),I(:,:,i)*w(:,i));
end
        
        

% 内推
f=zeros(3,n+1); %最后一列为外力，工具坐标系
t=zeros(3,n+1); %最后一列为外力矩，工具坐标系
for i=n:-1:1
    f(:,i)=R(:,:,i+1)*f(:,i+1)+F(:,i);
    t(:,i)=N(:,i)+R(:,:,i+1)*t(:,i+1)+cross(Pc(:,i),F(:,i))+cross(P(:,i+1),R(:,:,i+1)*f(:,i+1));
    if Links(i).sigma == 0
        Tau(i)=t(:,i)'*Z;
    else
        Tau(i)=f(:,i)'*Z;
    end        
end
end

