step=0.01;
fid3 = fopen('C:\Users\hqi\Desktop\motion_w_Newton_9\motion_w\pt_v.txt');
comma  = char(' ');
pt_v = fscanf(fid3, ['%f', comma]);
fclose(fid3);
[m,n]=size(pt_v);
pt_v(m)
plot([0:step:m*step-step],pt_v);title('pt_ V');