data=dlmread('C:\Users\hqi\Desktop\RPY_lin_blend_lin_joint_cmd.txt','',15);
[m,n]=size(data);
subplot(311);
plot(1:m,data);
subplot(312);
plot(1:m,data);
subplot(313);
plot(1:m,data);