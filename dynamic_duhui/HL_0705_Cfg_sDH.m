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
% DESCRIPTION:sDH方法建立六轴Helios-0705机器人模型，同时也包括机器人的驱动信息，其中的参数
%               已经在函数中写好。
%
% AUTHOR : Hui Du,09/2016
%
% ABSTRACT：这是用sDH方法建立六轴Helios-0705机器人的数学模型,其中参数
% 
% INPUT： none
%
% OUTPUT: config            包涵机器人所有参数的结构体   
%% 运动学和动力学建模参数
function config = HL_0705_Cfg_sDH ()
deg = pi/180;
config.type=6;                                                     % 定义机器人类型，为六轴机器人
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 运动学参数
% d        连杆偏置
% a        连杆长度
% alpha    连杆扭角
% theta    连杆转角
% sigma    0 for revolut, 1 for prismatic
% mDH      0 for sDH, 1 for mDH
% offset   关节转角偏置
% qlim     关节位置范围[min,max]
% 
% 动力学参数
% m        连杆质量
% r        连杆质心在连杆坐标系中的位置1x3
% I        惯性张量，即连杆对于坐标原点在质心，各轴方向与连杆坐标系方向相同的坐标系的惯性张量
% Jm       电机转子转动惯量  可删
file='HL_0705.xlsx';%不同表格配置不同地址
MotorData = '电机减速机参数';
MechData = '结构参数';
StructData = xlsread(file,MechData,'B4:S9');
for i = 1:6
    Links(i).d = StructData(i,3)*1e-3;
    Links(i).a = StructData(i,2)*1e-3;
    Links(i).alpha = StructData(i,1)*deg;
    Links(i).offset = StructData(i,4)*deg;
    Links(i).theta = NaN;
    Links(i).sigma = StructData(i,5);
    Links(i).mDH = StructData(i,6);
    Links(i).m = StructData(i,7);
    Links(i).r = StructData(i,8:10)*1e-3;
    I1 = diag(StructData(i,11:13));
    I2 = [0,StructData(i,14:15);0,0,StructData(i,16);0,0,0];
    Links(i).I = (I1+I2+I2')*1e-6;
    Links(i).qlim = StructData(i,17:18)*deg;
end
config.Links = Links;
config.mDH = 0;                         % 0 for standard DH ,1 for modified DH
config.Struct=[1 1 1];                  %机器人姿态，和求解反解相关;1*3向量，取值只有1或者2
                                        % Struct(1)=1 front 正面
                                        % Struct(1)=2 back  反面
                                        % Struct(2)=1 up    肘部向上
                                        % Struct(2)=2 down  肘部向下
                                        % Struct(3)=1 no flip    腕关节J5正转
                                        % Struct(3)=2 fliped over    腕关节J5反转   
config.Load = Load_0705_Cfg(5);         % 配置5kg负载
config.Jm = xlsread(file,MotorData,'B20:G20')*1e-6;
%% 电机减速机相关运动参数
% 关节到电机的减速比
config.G = xlsread(file,MotorData,'B15:G15');
% 最大关节转速
config.VmaxMat = xlsread(file,MotorData,'B6:G6')./config.G*6*pi/180; 
config.AmaxMat = config.VmaxMat*10;                                % 最大关节角加速度

config.VmaxCartesian=[7 4*pi];                                     % 最大笛卡尔空间速度，[合成线速度，合成角速度]
config.AmaxCartesian=config.VmaxCartesian*10;                      % 最大笛卡尔空间加速度，[合成线加速度，合成角加速度]

 % 关节转速极限
config.VLimit = xlsread(file,MotorData,'B18:G18')./config.G*6*pi/180;
config.ALimit=config.VLimit* 10;                                   % 关节角加速度极限

% 关节力矩极限
config.TLimit = xlsread(file,MotorData,'B19:G19');
config.MotorTorque = xlsread(file,MotorData,'B5:G5');

config.Kt = xlsread(file,MotorData,'B9:G9');                         % 电机力矩系数，单位N*m/A
config.R = xlsread(file,MotorData,'B10:G10');                        % 电机相电阻   

config.eff = 0.9;                                                    % 功率转换系数
config.MechEff = ones(1,6);                                          % 机械效率
config.Tf = [0.00127,0.00127,0.000667,0.000667,0.000267 ,0.000267];  % 粘滞摩擦系数   

%% 仿真参数
config.StepSize = 0.001;                                             % 仿真步长    

