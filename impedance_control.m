clear;
l1=0.25;l2=0.25;  m1=2.247;m2=11.078;  %两连杆长度质量
g=0;%重力加速度
Md=[1 0;0 1]; Kd=[2 0;0 2]; Kp=[1 0;0 1];

t=[1 ;2];%力矩
te=[10;10];%外力
theta=[10 ;10]; 

theta_d=[0 ;0];
ddtheta_d=[0;0];
dtheta=[0; 0];
ddtheta=[0 ; 0];
times=1000;
plot_theta=zeros(times,2);
plot_t=zeros(times,2);
for i=1:times
M=[l2^2*m2 + 2*l1*l2*m2*cos(theta(2)) + l1^2*(m1+m2)  l2^2*m2 + l1*l2*m2*cos(theta(2))
    l2^2*m2 + l1*l2*m2*cos(theta(2))                  l2^2*m2];
C=[-m2*l1*l2*sin(theta(2))*dtheta(2)^2 - 2*m2*l1*l2*sin(theta(2))*dtheta(1)*dtheta(2)
    m2*l1*l2*sin(theta(2))*dtheta(1)^2];
G=[m2*l2*g*cos(theta(1)+theta(2)) + (m1+m2)*l1*g*cos(theta(1))
    m2*l2*g*cos(theta(1)+theta(2))];
F=0.2*[1 0;0 2]*dtheta+[0.5 ; 0.6].*(abs(dtheta)>1E-10)+[1 ; 1].*(abs(dtheta)<=1E-10);
%J=[-l1*sin(theta(1))-l2*sin(theta(1)+theta(2))  -l2*sin(theta(1)+theta(2))
%    l1*cos(theta(1))+l2*cos(theta(1)+theta(2))  l2*cos(theta(1)+theta(2))];%雅克比矩阵
t=M*(0+8*[1 0;0 1]*(0-dtheta)+16*(theta_d-theta))+C+G+F;%-->>PD
%t=t+te; %-->>PD
% y=9*Md^-1*(Kd*-dtheta+Kp*-theta)-dtheta;
% t=M*y+C+G;
plot_t(i,:)=t;
ddtheta=M^-1*(t-C-G-F);
dtheta=dtheta+ddtheta.*0.01;
theta=theta+dtheta.*0.01+0.5*ddtheta.*(0.01)^2;
plot_theta(i,:)=theta(:);

end
subplot(211);plot(1:times,plot_theta(:,1),1:times,plot_theta(:,2));
subplot(212);plot(1:times,plot_t);