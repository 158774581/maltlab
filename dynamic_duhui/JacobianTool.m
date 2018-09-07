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
% ABSTRACT����⽫�ؽ��ٶ�ӳ��Ϊĩ��ִ�����ڹ�������ϵ�е��ٶȵ������ſ˱Ⱦ���ͨ�ú���
% 
% INPUT��q            �����˹ؽ�λ�Ʊ�������λrad��m��1xN��Nx1����
%        config       ���������ýṹ��
%        
% OUTPUT: J           �����ſ˱Ⱦ���   
%                     config.type = 4ʱ��scara�����ˣ�JΪ4x4����
%                     config.type = 6ʱ����������ˣ�JΪ6x6����
% 
function J = JacobianTool(q,config)
L = config.Links;
U = eye(4); %%������
n = config.type;
for i=n:-1:1
    if config.mDH == 0
        % standard DH
        U = Transformation(q(i),L(i)) * U;
    end
    if L(i).sigma == 0
        % revolute joint
        d = [-U(1,1)*U(2,4)+U(2,1)*U(1,4)
             -U(1,2)*U(2,4)+U(2,2)*U(1,4)
             -U(1,3)*U(2,4)+U(2,3)*U(1,4)];        
        delta = U(3,1:3)';
    else
        % prismatic joint 
        d = U(3,1:3)';
        delta = [0;0;0];
    end
    J(:,i) = [d;delta];
    if config.mDH == 1
        U = Transformation(q(i),L(i)) * U;
    end
end

if n == 4
    % for SCARA
    J = [J(1:3,:);J(6,:)];
end
        
        
        
        
