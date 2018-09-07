% step=0.001;
% fid4 = fopen('C:\Users\hqi\Desktop\llinear_equation\motion_w\abcd.txt');
% comma  = char(' ');
% abcd = fscanf(fid4, ['%f', comma]);
% fclose(fid4);
% [m,n]=size(abcd);
A=[
1.000000  1.000000  1.000000  1.000000  0.475456  
3.000000  2.000000  1.000000  0.000000  1.247870  
1.003003  1.002001  1.001000  1.000000  0.476705 
3.006003  2.002000  1.000000  0.000000  1.249139 
 ];
[L,u]=lu(A);
b=[4 4 4 4];
x=A/b
% A=[1 1 1 2 3;3 2 1 4 5 ; 1 2 2 3 2; 4 2 4 1 4];
% [L,v]=lu(A);