% inertia
% ����SCARA������1,2���ת������.

% ���ֹ�ʽ:
% ����J1Ϊ�����ۺϵ�����Ĳ�ת�������͵������ת�������ĺ�.
% ��ΪJ2����̬��ͬ������J1��ת��������ͬ,������ͬ��λ����һ��ƽ��ֵ,��Ϊ1�ؽڼ��ٱ�reduction_ratio_1
% J1_load_to_motor= (J1_inertia_152+J1_inertia_90+J1_inertia_0)/3/reduction_ratio_1^2;
% J1=J1_load_to_motor+J1_motor
% Ϊ�˷�ֹ��ת�����Ƚ�С��ת�������ĵط��ٶȻ�����̫����ʱѡȡ�Ƚ�С�Ĺ���ֵ��ƽ��ֵ��2/3
% 
% J2�ļ������ơ�

% INPUTS : m                    ��������     
%          reduction_ratio_1    J1���ٱ�
%          reduction_ratio_2    J1���ٱ�
%          J1_motor             J1�������ת������
%          J2_motor             J2�������ת������
% OUTPUTS :
%           J1                  1���ת������
%           J1_067              1���ת��������2/3
%           J2                  2���ת������


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
        %��J2��ת�����������ֵ
        J2_inertia=0.00259*m^3 - 0.02044*m^2 +  0.08714*m + 0.2063;  
end
J2_load_to_motor=J2_inertia/reduction_ratio_2^2;
J2=J2_load_to_motor+J2_motor;
disp('J2ת������')
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
        %��J2��ת�����������ֵ
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
        %��J2��ת�����������ֵ
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
        %��J2��ת�����������ֵ
        J1_inertia_152= 0.0003906*m^3-0.002992*m^2 + 0.01513*m +  0.368;
end
J1_load_to_motor =(J1_inertia_152+J1_inertia_90+J1_inertia_0)/3/reduction_ratio_1^2;
J1=J1_load_to_motor+J1_motor;
disp('J1ת������')
disp(J1);
J1_067=2/3*J1;
disp('J1ת��������2/3')
disp(J1_067);

end
%��������
%m_load=[0 2 5 6];
%��J2��ת������
%inertia_J2=[0.206317 0.319568 0.454866 0.552914];
%��J2��ת�����������ֵ
%inertia_J2_p=0.00259*m_load.^3 - 0.02044*m_load.^2 +  0.08714*m_load + 0.2063;





