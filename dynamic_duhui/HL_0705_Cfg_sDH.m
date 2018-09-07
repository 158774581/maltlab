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
% DESCRIPTION:sDH������������Helios-0705������ģ�ͣ�ͬʱҲ���������˵�������Ϣ�����еĲ���
%               �Ѿ��ں�����д�á�
%
% AUTHOR : Hui Du,09/2016
%
% ABSTRACT��������sDH������������Helios-0705�����˵���ѧģ��,���в���
% 
% INPUT�� none
%
% OUTPUT: config            �������������в����Ľṹ��   
%% �˶�ѧ�Ͷ���ѧ��ģ����
function config = HL_0705_Cfg_sDH ()
deg = pi/180;
config.type=6;                                                     % ������������ͣ�Ϊ���������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �˶�ѧ����
% d        ����ƫ��
% a        ���˳���
% alpha    ����Ť��
% theta    ����ת��
% sigma    0 for revolut, 1 for prismatic
% mDH      0 for sDH, 1 for mDH
% offset   �ؽ�ת��ƫ��
% qlim     �ؽ�λ�÷�Χ[min,max]
% 
% ����ѧ����
% m        ��������
% r        ������������������ϵ�е�λ��1x3
% I        ���������������˶�������ԭ�������ģ����᷽������������ϵ������ͬ������ϵ�Ĺ�������
% Jm       ���ת��ת������  ��ɾ
file='HL_0705.xlsx';%��ͬ������ò�ͬ��ַ
MotorData = '������ٻ�����';
MechData = '�ṹ����';
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
config.Struct=[1 1 1];                  %��������̬������ⷴ�����;1*3������ȡֵֻ��1����2
                                        % Struct(1)=1 front ����
                                        % Struct(1)=2 back  ����
                                        % Struct(2)=1 up    �ⲿ����
                                        % Struct(2)=2 down  �ⲿ����
                                        % Struct(3)=1 no flip    ��ؽ�J5��ת
                                        % Struct(3)=2 fliped over    ��ؽ�J5��ת   
config.Load = Load_0705_Cfg(5);         % ����5kg����
config.Jm = xlsread(file,MotorData,'B20:G20')*1e-6;
%% ������ٻ�����˶�����
% �ؽڵ�����ļ��ٱ�
config.G = xlsread(file,MotorData,'B15:G15');
% ���ؽ�ת��
config.VmaxMat = xlsread(file,MotorData,'B6:G6')./config.G*6*pi/180; 
config.AmaxMat = config.VmaxMat*10;                                % ���ؽڽǼ��ٶ�

config.VmaxCartesian=[7 4*pi];                                     % ���ѿ����ռ��ٶȣ�[�ϳ����ٶȣ��ϳɽ��ٶ�]
config.AmaxCartesian=config.VmaxCartesian*10;                      % ���ѿ����ռ���ٶȣ�[�ϳ��߼��ٶȣ��ϳɽǼ��ٶ�]

 % �ؽ�ת�ټ���
config.VLimit = xlsread(file,MotorData,'B18:G18')./config.G*6*pi/180;
config.ALimit=config.VLimit* 10;                                   % �ؽڽǼ��ٶȼ���

% �ؽ����ؼ���
config.TLimit = xlsread(file,MotorData,'B19:G19');
config.MotorTorque = xlsread(file,MotorData,'B5:G5');

config.Kt = xlsread(file,MotorData,'B9:G9');                         % �������ϵ������λN*m/A
config.R = xlsread(file,MotorData,'B10:G10');                        % ��������   

config.eff = 0.9;                                                    % ����ת��ϵ��
config.MechEff = ones(1,6);                                          % ��еЧ��
config.Tf = [0.00127,0.00127,0.000667,0.000667,0.000267 ,0.000267];  % ճ��Ħ��ϵ��   

%% �������
config.StepSize = 0.001;                                             % ���沽��    

