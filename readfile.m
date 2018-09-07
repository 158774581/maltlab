fid = fopen('C:\Users\hqi\Documents\Visual Studio 2017\Projects\motion_PVT_linear\motion_w_Newton_9\motion_w\a.txt');
comma  = char(' ');
a = fscanf(fid, ['%f', comma]);
fclose(fid);
[m,n]=size(a);
disp('末加速度'); 
a(m)
step=0.001;
fid1 = fopen('C:\Users\hqi\Documents\Visual Studio 2017\Projects\motion_PVT_linear\motion_w_Newton_9\motion_w\v.txt');
comma  = char(' ');
v = fscanf(fid1, ['%f', comma]);
fclose(fid1);
disp('末速度');
v(m)
fid2 = fopen('C:\Users\hqi\Documents\Visual Studio 2017\Projects\motion_PVT_linear\motion_w_Newton_9\motion_w\d.txt');
comma  = char(' ');
d = fscanf(fid2, ['%f', comma]);
disp('末位移');
d(m)
fclose(fid2);
subplot(311);
plot([0:step:m*step-step],a);title('加速度');
subplot(312);
plot([0:step:m*step-step],v);title('速度');
subplot(313);
plot([0:step:m*step-step],d);title('位移');