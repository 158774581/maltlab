% 
% syms m1 t m2 l1 l2 h1 h2 w1 w2 p1 p2 v1 v2 a1 a2 c1x c1y c1z c2x c2y c2z real;
% syms p_1_c1 p_2_c2 I_c1_1 I_c2_2  real;
% p_1_c1=[c1x;0;0];%input
% p_2_c2=[c2x;0;0];%input
% syms I1xx I1xy I1xz I1yy I1yz I1zz  I2xx I2xy I2xz I2yy I2yz I2zz real;
% I_c1_1=[I1xx -0 -0;-0 I1yy -0;-0 -0 I1zz]; %input
% I_c2_2=[I2xx -0 -0;-0 I2yy -0;-0 -0 I2zz]; %input
% w_0=0; w_d_0=0; f_3=0; n_3=0; V_d_0_0=0;%gravity
% R_0_1=[cos(p1) -sin(p1) 0;sin(p1) cos(p1) 0;0 0 1];
% R_1_0=R_0_1';
% R_1_2=[cos(p2) -sin(p2) 0;sin(p2) cos(p2) 0;0 0 1];
% R_2_1=R_1_2';
% %------------开始外推 1 ---------------
% w_1_1=v1*[0;0;1];  w_1_1_m=[0 -v1 0;v1 0 0;0 0 0];
% w_d_1_1=[0;0;a1];
% v_d_1_1=0;
% v_d_1_c1=cross(w_d_1_1,p_1_c1)+cross(w_1_1,cross(w_1_1,p_1_c1))+v_d_1_1;
% F_1_1=m1*v_d_1_c1;
% N_1_1=I_c1_1*w_d_1_1+cross(w_1_1,I_c1_1*w_1_1);
% %------------开始外推 2 ---------------
% w_2_2=R_2_1*w_1_1+v2*[0;0;1];
% w_2_2_m=[0 -w_2_2(3) w_2_2(2);w_2_2(3) 0 -w_2_2(1);-w_2_2(2) w_2_2(1) 0];
% w_d_2_2=R_2_1*w_d_1_1+R_2_1*cross(w_1_1,[0;0;v2])+[0;0;a2];
% v_d_2_2=R_2_1*(cross(w_d_1_1,[ l1;0;0])+cross(w_1_1,cross(w_1_1,[l1;0;0]))+v_d_1_1);
% v_d_2_c2=cross(w_d_2_2,p_2_c2)+cross(w_2_2,cross(w_2_2,p_2_c2))+v_d_2_2;
% F_2_2=m2*v_d_2_c2;
% N_2_2=I_c2_2*w_d_2_2+cross(w_2_2,I_c2_2*w_2_2);
% %------------开始内推 2 ---------------
% %                                                                                                                                                                                                                                                                                                                                                   
% f_2_2=0+F_2_2;
% n_2_2=N_2_2+0+cross(p_2_c2,F_2_2);
% t2=n_2_2'*[0;0;1];
% %------------开始内推 1 ---------------
% P_1_2=[l1;0;0];
% P_1_2_m=[0 -P_1_2(3) P_1_2(2);P_1_2(3) 0 -P_1_2(1);-P_1_2(2) P_1_2(1) 0];
% f_1_1=R_1_2*f_2_2+F_1_1;    
% n_1_1=N_1_1+R_1_2*n_2_2+cross(p_1_c1,F_1_1)+cross(P_1_2,R_1_2*f_2_2);
% t1=n_1_1'*[0;0;1];
% t1=collect(t1,p1)
% t1=collect(t1,v1);
% t1=collect(t1,a1);
% t2=collect(t2,p2);
% t2=collect(t2,v2) 
% t2=collect(t2,a2);
% 
% 
joint_fd=dlmread('C:\Users\hqi\Desktop\impedance_joint.txt','',0);
joint_d=dlmread('C:\Users\hqi\Desktop\joint_cmd.txt','',0);
cart_fd=dlmread('C:\Users\hqi\Desktop\impedance_cart.txt','',0);
cart_d=dlmread('C:\Users\hqi\Desktop\cart_cmd.txt','',0);
torque=dlmread('C:\Users\hqi\Desktop\impedance_torque.txt','',0);
aq=dlmread('C:\Users\hqi\Desktop\备用.txt','',0);
[m,n]=size(joint_fd);
[m1,n1]=size(joint_d);
joint_d=[joint_d;zeros(m-m1,6)];
cart_d=[cart_d;zeros(m-m1,6)];
subplot(211);
plot(1:m,joint_fd(:,1),1:m,joint_d(:,1),1:m,joint_fd(:,2),1:m,joint_d(:,2));
subplot(212);
plot(1:m,cart_fd(:,1),1:m,cart_d(:,1),1:m,cart_fd(:,2),1:m,cart_d(:,2));
figure(2);
subplot(211);
plot(1:m,torque(1:m,1),'r',1:m,torque(1:m,2),'g');
subplot(212);
plot(1:m,torque(1:m,3),'r',1:m,torque(1:m,4),'g');
figure(3);
plot(aq);

% flag=0;
% while(flag==0)
% r=unifrnd(-5,5,10,1);
% P=[r(1) r(2)  r(3) r(4) 
%   r(2)  r(5)  r(6) r(7)
%   r(3) r(6)  r(8)  r(9)
%   r(4) r(7)  r(9)  r(10)];
% 
% 
% 
% A=[ 0     0     1     0
%      0     0     0     1
%     -1     0    -1     0
%      0    -1     0    -1];
%  Q=-P*A-A'*P;
%  Q_eig=eig(Q);
%  P_eig=eig(P);
%  flag=1;
%  for i=1:4
%  if(Q_eig(i)<0||isreal(Q_eig(i))~=1||P_eig(i)<0||isreal(P_eig(i))~=1)
%  flag=0;
%  end
%  end
% end
%  P
%  Q
%  P_eig
%  Q_eig


 