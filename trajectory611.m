clear all;
joint_cmd=dlmread('C:\Users\hqi\Desktop\joint_cmd.txt','',0);
cart_cmd=dlmread('C:\Users\hqi\Desktop\cart_cmd.txt','',0);
cart_vel_cmd=dlmread('C:\Users\hqi\Desktop\cart_vel_cmd.txt','',0);
joint_vel_cmd=dlmread('C:\Users\hqi\Desktop\joint_vel_cmd.txt','',0);
 
[m,n]=size(joint_cmd);
[l,p]=size(cart_cmd);
[o,q]=size(cart_vel_cmd);
[m1.n1]=size(joint_vel_cmd);

figure(1);
subplot(311); %joint pos
plot(1:m,joint_cmd(:,1));title('joint-cmd');
hold on;
plot(1:m,joint_cmd(:,2));
hold on;
plot(1:m,joint_cmd(:,3));
hold on;
plot(1:m,joint_cmd(:,4));
hold on;
plot(1:m,joint_cmd(:,5));
hold on;
plot(1:m,joint_cmd(:,6));
hold on;
legend('j1','j2','j3','j4','j5','j6');
hold off;

subplot(312); %cart pos
plot(1:l,cart_cmd(:,1));title('cart12-cmd');
hold on;
plot(1:l,cart_cmd(:,2));
hold on;
subplot(313); %cart pos
 plot(1:l,cart_cmd(:,3));title('cart-3456-cmd');
hold on;   
plot(1:l,cart_cmd(:,4));
hold on;
plot(1:l,cart_cmd(:,5));
hold on;
plot(1:l,cart_cmd(:,6));
hold on;
legend('c1','c2','c3','c4','c5','c6');
hold off;
figure(3)
subplot(211);plot(1:o,joint_vel_cmd(:,1:2));  title('joint12-vel');
subplot(212);plot(1:o,joint_vel_cmd(:,3:6)); title('joint3456-vel');
joint=zeros(6,1);%joint空间的位置
for i=1:6
 joint(i)=sum(joint_vel_cmd(:,i))*0.004;
end
joint
figure(4);%cart vel
subplot(211);plot(1:o,cart_vel_cmd(:,1:2));  title('cart12-vel');
subplot(212);plot(1:o,cart_vel_cmd(:,3:6)); title('cart3456-vel');
cart=zeros(6,1); %笛卡尔坐标下的位移、姿态
for i=1:6
 cart(i)=sum(cart_vel_cmd(:,i))*0.004;
end
cart