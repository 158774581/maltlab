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
% ABSTRACT����⽫�ؽ��ٶ�ӳ��Ϊĩ��ִ�����ڿռ�����ϵ�е��ٶȵĿռ��ſ˱Ⱦ���ͨ�ú���
%           ��Jacobian����������ͬ��������ͬ
% 
% INPUT��q            �����˹ؽ�λ�Ʊ�������λrad��m��1xN��Nx1����
%        config       ���������ýṹ��
%        
% OUTPUT: J           �ռ��ſ˱Ⱦ���   
%                     config.type = 4ʱ��scara�����ˣ�JΪ4x4����
%                     config.type = 6ʱ����������ˣ�JΪ6x6����
% 
function J = JacobianWorld(q,config)
Jn = JacobianTool(q,config);
Tn = Fkine(q,config);
R = t2r(Tn);
n = config.type;
if n == 6 || n == 7
    % ����
    J = [R zeros(3);zeros(3) R]*Jn;
elseif n == 4
    % for SCARA 
    J = [R zeros(3,1);0,0,0,R(3,3)]*Jn;
end
    











