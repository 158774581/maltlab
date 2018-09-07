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
% ABSTRACT：求解将关节速度映射为末端执行器在空间坐标系中的速度的空间雅克比矩阵，通用函数
%           与Jacobian函数功能相同，方法不同
% 
% INPUT：q            机器人关节位移变量，单位rad或m，1xN或Nx1向量
%        config       机器人配置结构体
%        
% OUTPUT: J           空间雅克比矩阵   
%                     config.type = 4时，scara机器人，J为4x4矩阵
%                     config.type = 6时，六轴机器人，J为6x6矩阵
% 
function J = JacobianWorld(q,config)
Jn = JacobianTool(q,config);
Tn = Fkine(q,config);
R = t2r(Tn);
n = config.type;
if n == 6 || n == 7
    % 六轴
    J = [R zeros(3);zeros(3) R]*Jn;
elseif n == 4
    % for SCARA 
    J = [R zeros(3,1);0,0,0,R(3,3)]*Jn;
end
    











