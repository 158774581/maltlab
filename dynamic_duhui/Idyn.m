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
% DESCRIPTION:��ţ��-ŷ����������ѧ�㷨���ɹؽ��˶�����ؽ����ء���һ�����Ƕ�ÿ��
%             ����Ӧ��ţ��-ŷ�����̣�������1������n��������������˵��ٶȺͼ��ٶȡ�
%             �ڶ������Ǵ�����n������1���ڵ����������˼���໥�������������Լ��ؽ�
%             �������ء�
%
% AUTHOR : Hui Du,09/2016
%
% ABSTRACT������Helios����������涯��ѧ������ͨ�ú����������ڴ���������
%           ��ţ��-ŷ����������ѧ�㷨���ɹؽ��˶�����ؽ����ء���һ�����Ƕ�ÿ��
%           ����Ӧ��ţ��-ŷ�����̣�������1������n��������������˵��ٶȺͼ��ٶȡ�
%           �ڶ������Ǵ�����n������1���ڵ����������˼���໥�������������Լ��ؽ�
%           �������ء�
% 
% INPUT��q                �ؽڽǶȣ�1xN���󣬵�λrad
%        qd               �ؽڽ��ٶȣ�1xN����λrad/s
%        qdd              �ؽڽǼ��ٶȣ�1xN���󣬵�λrad/(s^2)
%        config           ���������ýṹ��
%
% OUTPUT: Tau             �ؽ����ؾ���1xN����λNm       
% 
% REFERENCE����������ѧ���ۡ� 
%
function Tau = Idyn( q,qd,qdd,config)
% setting 
Links=config.Links;
n=config.type;
% 1-n�����˵���������
for i=1:n-1
     m(i)=Links(i).m;
    if Links(i).sigma == 1
        Links(i).d = q(i);
        q(i) = 0; %%%����֤��һ��
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

%ĩ��ִ�������������ԣ��Ӹ���
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
% ����
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
        
        

% ����
f=zeros(3,n+1); %���һ��Ϊ��������������ϵ
t=zeros(3,n+1); %���һ��Ϊ�����أ���������ϵ
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

