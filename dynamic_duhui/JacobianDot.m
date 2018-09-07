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
% ABSTRACT������ſ˱Ⱦ������ĺ�����ͨ�ú���
%           ��������JdotΪ��Ӧ��ȫ������ϵ�µ�
% 
% INPUT��q            �����˹ؽ�λ�Ʊ�������λrad��m��1xN��Nx1����
%        qd           �����˹ؽ��ٶȱ�������λrad/s��m/s��1xN��Nx1����
%        config       ���������ýṹ��
%        
% OUTPUT: Jdot        �ſ˱Ⱦ���ĵ���  
%                     config.type = 4ʱ��scara�����ˣ�JΪ4x4����
%                     config.type = 6ʱ����������ˣ�JΪ6x6����
%
function Jdot = JacobianDot(q,qd,config)

%   T       ÿ�����˵Ĵ��ݺ���/����任����4x4����
%   Q��a    ÿ�����˵���ת���󻻾����λ��������3x3�����3x1����
%   P       ȫ������ϵ����������ϵ����ת����3x3����
%   e       ȫ������ϵ�У��ؽ����ߵķ���������3x1����
%   w       ���˽��ٶȣ�3x1����
%   ed      e��΢�֣�3x1����
%   ri      ��ǰ����ϵԭ�㵽ĩ�˹�������ϵԭ����������ڵ�ǰ��������ϵ�б�ʾ��3x1����
%   r0      ri��ȫ������ϵ�б�ʾ��3x1����
%   a0      a��ȫ������ϵ�еı�ʾ���꣬3x1����
%   rd      r0��΢�֣�3x1����
%   Bd      �ſ˱Ⱦ������ƶ�������΢�֣�3x6����
%   ed      �ſ˱Ⱦ�����ת��������΢�֣�3x6����


n = config.type;
L = config.Links;
Q = zeros(3,3,n);
P = zeros(3,3,n);
e = zeros(3,n);
a = zeros(3,n);

% �������˴��ݺ���
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

 








    
    



    