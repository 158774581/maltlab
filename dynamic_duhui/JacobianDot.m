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
% AUTHOR : Hui Du, 09/2016
%
% ABSTRACT：求解雅克比矩阵导数的函数，通用函数
%           但求解出的Jdot为对应到全局坐标系下的
% 
% INPUT：q            机器人关节位移变量，单位rad或m，1xN或Nx1向量
%        qd           机器人关节速度变量，单位rad/s或m/s，1xN或Nx1向量
%        config       机器人配置结构体
%        
% OUTPUT: Jdot        雅克比矩阵的导数  
%                     config.type = 4时，scara机器人，J为4x4矩阵
%                     config.type = 6时，六轴机器人，J为6x6矩阵
%
function Jdot = JacobianDot(q,qd,config)

%   T       每个连杆的传递函数/坐标变换矩阵，4x4矩阵
%   Q，a    每个连杆的旋转矩阵换矩阵和位置向量，3x3矩阵和3x1向量
%   P       全局坐标系到连杆坐标系的旋转矩阵，3x3矩阵
%   e       全局坐标系中，关节轴线的方向向量，3x1向量
%   w       连杆角速度，3x1向量
%   ed      e的微分，3x1向量
%   ri      当前坐标系原点到末端工具坐标系原点的向量，在当前连杆坐标系中表示，3x1向量
%   r0      ri在全局坐标系中表示，3x1向量
%   a0      a在全局坐标系中的表示坐标，3x1向量
%   rd      r0的微分，3x1向量
%   Bd      雅克比矩阵中移动分量的微分，3x6向量
%   ed      雅克比矩阵中转动分量的微分，3x6向量


n = config.type;
L = config.Links;
Q = zeros(3,3,n);
P = zeros(3,3,n);
e = zeros(3,n);
a = zeros(3,n);

% 计算连杆传递函数
for i = 1:n
    T = Transformation(q(i),L(i));
    Q(:,:,i) = t2r(T);
    a(:,i) = T(1:3,4);
end

P(:,:,1) = Q(:,:,1);
e(:,1) = [0;0;1];

if config.mDH == 0
    for i = 2:n
        P(:,:,i) = P(:,:,i-1) * Q(:,:,i);
        e(:,i) =  P(:,3,i-1);
    end
    
    % step 1
    w(:,1) = qd(1) * e(:,1);
    for i = 1:n-1
        if L(i+1).sigma == 0
            w(:,i+1) = qd(i+1)*e(:,i+1)+w(:,i); %world frame
        else
            w(:,i+1) = w(:,i);
        end
    end
    
    % step 2
    ed(:,1) = zeros(3,1);
    for i = 2:n
        ed(:,i) = cross(w(:,i),e(:,i));
    end
    
    % step 3
    ri(:,n) = a(:,n);
    for i = n-1:-1:1
        ri(:,i) = a(:,i) + Q(:,:,i)*ri(:,i+1);
    end
    
    r0(:,1) = ri(:,1);
    for i = 2:n
        r0(:,i) = P(:,:,i-1)*ri(:,i);
    end
    
    a0(:,1) = a(:,1);
    for i = 2:n
        a0(:,i) = P(:,:,i-1)*a(:,i);
    end
    
    rd(:,n) = cross(w(:,n),a0(:,n));
    for i = n-1:-1:1
        rd(:,i) = cross(w(:,i),a0(:,i))+rd(:,i+1);
    end
    
    % step 4
    Bd(:,1) = cross(e(:,1),rd(:,1));
    for i = 2:n
        Bd(:,i) = cross(ed(:,i),r0(:,i))+cross(e(:,i),rd(:,i));
    end

else
    for i = 2:n
        P(:,:,i) = P(:,:,i-1) * Q(:,:,i);
        e(:,i) =  P(:,3,i);
    end
    Q(:,:,n+1) = eye(3);
    
    % step 1
    w(:,1) = qd(1) * e(:,1);
    for i = 1:n-1
        if L(i+1).sigma == 0
            w(:,i+1) = qd(i+1)*e(:,i+1)+w(:,i); %world frame
        else
            w(:,i+1) = w(:,i);
        end
    end
    
    % step 2
    ed(:,1) = zeros(3,1);
    for i = 2:n
        ed(:,i) = cross(w(:,i),e(:,i));
    end
    
    % step 3
    ri(:,n) = [0;0;0];
    for i = n-1:-1:1
        ri(:,i) = a(:,i+1) + Q(:,:,i+1)*ri(:,i+1);
    end
    ri = [a(:,1) + Q(:,:,1)*ri(:,1),ri];
    
    r0(:,1) = ri(:,1);
    for i = 2:n+1
        r0(:,i) = P(:,:,i-1)*ri(:,i);
    end
    
    a0(:,1) = a(:,1);
    a(:,n) = [0;0;0];
    for i = 2:n+1
        a0(:,i) = P(:,:,i-1)*a(:,i); 
    end
    
    rd(:,n) = cross(w(:,n),a0(:,n+1));
    for i = n-1:-1:1
        rd(:,i) = cross(w(:,i),a0(:,i+1))+rd(:,i+1);
    end
    
    % step 4
    Bd(:,1) = cross(e(:,1),rd(:,1));
    for i = 2:n
        Bd(:,i) = cross(ed(:,i),r0(:,i+1))+cross(e(:,i),rd(:,i));
    end
        
end
% step 5
Jdot = [];
for i = 1:n
    if L(i).sigma == 0
        Jdot = [Jdot,[Bd(:,i);ed(:,i)]];
    else
        Jdot = [Jdot,[ed(:,i);0;0;0]];
    end
end

if n == 4
    Jdot = [Jdot(1:3,:);Jdot(6,:)];
end

 








    
    



    