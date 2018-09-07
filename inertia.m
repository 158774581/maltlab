% inertia
% 计算SCARA机器人1,2轴的转动惯量.

% 部分公式:
% 整个J1为负载折合到电机的侧转动惯量和电机本身转动惯量的和.
% 因为J2的姿态不同，导致J1的转动惯量不同,几个不同的位置求一个平均值,因为1关节减速比reduction_ratio_1
% J1_load_to_motor= (J1_inertia_152+J1_inertia_90+J1_inertia_0)/3/reduction_ratio_1^2;
% J1=J1_load_to_motor+J1_motor
% 为了防止在转动到比较小的转动惯量的地方速度环增益太大，暂时选取比较小的惯量值，平均值的2/3
% 
% J2的计算类似。

% INPUTS : m                    负载质量     
%          reduction_ratio_1    J1减速比
%          reduction_ratio_2    J1减速比
%          J1_motor             J1电机本身转动惯量
%          J2_motor             J2电机本身转动惯量
% OUTPUTS :
%           J1                  1轴的转动惯量
%           J1_067              1轴的转动惯量的2/3
%           J2                  2轴的转动惯量


function [J1,J1_067,J2]=inertia(m,reduction_ratio_1,reduction_ratio_2,J1_motor,J2_motor)
format long ;
switch m
    case 0
        J2_inertia=0.206317;
    case 2
        J2_inertia=0.319568;
    case 5
        J2_inertia=0.454866;
    case 6
        J2_inertia=0.552914;
    otherwise
        %绕J2的转动惯量的拟合值
        J2_inertia=0.00259*m^3 - 0.02044*m^2 +  0.08714*m + 0.2063;  
end
J2_load_to_motor=J2_inertia/reduction_ratio_2^2;
J2=J2_load_to_motor+J2_motor;
disp('J2转动惯量')
disp(J2);
switch m
      case 0
        J1_inertia_0=1.384126;
    case 2
        J1_inertia_0=1.878867;
    case 5
        J1_inertia_0=2.2484885;
    case 6
        J1_inertia_0=2.863431;
    otherwise
        %绕J2的转动惯量的拟合值
        J1_inertia_0=0.02463*m^3 -0.1972*m^2 +  0.5433*m + 1.384;
end
switch m
      case 0
        J1_inertia_90=0.840651;
    case 2
        J1_inertia_90=1.084057;
    case 5
        J1_inertia_90=1.327436;
    case 6
        J1_inertia_90=1.570035;
    otherwise
        %绕J2的转动惯量的拟合值
        J1_inertia_90= 0.008081*m^3-0.06468*m^2 +  0.2187*m + 0.8407;
end
switch m
      case 0
        J1_inertia_152=0.368028;
    case 2
        J1_inertia_152=0.389448;
    case 5
        
        J1_inertia_152=0.417709;
    case 6
        J1_inertia_152=0.435472;
    otherwise
        %绕J2的转动惯量的拟合值
        J1_inertia_152= 0.0003906*m^3-0.002992*m^2 + 0.01513*m +  0.368;
end
J1_load_to_motor =(J1_inertia_152+J1_inertia_90+J1_inertia_0)/3/reduction_ratio_1^2;
J1=J1_load_to_motor+J1_motor;
disp('J1转动惯量')
disp(J1);
J1_067=2/3*J1;
disp('J1转动惯量的2/3')
disp(J1_067);

end
%负载质量
%m_load=[0 2 5 6];
%绕J2的转动惯量
%inertia_J2=[0.206317 0.319568 0.454866 0.552914];
%绕J2的转动惯量的拟合值
%inertia_J2_p=0.00259*m_load.^3 - 0.02044*m_load.^2 +  0.08714*m_load + 0.2063;





