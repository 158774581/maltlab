% %结合在一起的总程序
% %motioncase 1
% clear;%速度达到最大，有匀速运动
% figure(1);
% J1=3;J2=2;
% step=0.001;
% Va=6;   %机器允许的最大速度
% S=22;
% Acc=3;
% Dec=2;
% Tacc=ceil(Acc/(J1*step))*step;
% Tdec=ceil(Dec/(J2*step))*step;
% J1=Acc/Tacc;
% J2=Dec/Tdec;
% T1=ceil((Va/(Tacc*J1)-Tacc)/step)*step;
% J1=Va/(Tacc^2+T1*Tacc);
% T2=ceil((Va/(Tdec*J2)-Tdec)/step)*step;
% J2=Va/(Tdec^2+T2*Tdec);
% 
% C1=Va*Tacc+Va*T1/2+Va*Tdec+Va*T2/2;
% T3=ceil((S/(Tacc^2+T1*Tacc)/J1-(Tacc+T1/2+Tdec+T2/2))/step)*step;
% T=2*(Tacc+Tdec)+T1+T2+T3;
% 
% t=0:step:T;
% a=zeros(1,ceil(T/step)-1);
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% v=zeros(size(a));
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% d=zeros(size(a));
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
% 
%motioncase 2
% clear;%速度没有达到最大，无匀速运动。加（减）速度达到最大
% figure(2);
% J1=3;J2=2;
% step=0.01;
% Va=6;   %机器允许的最大速度
% S=18;
% Acc=3;
% Dec=2;
% Tacc=ceil(Acc/(J1*step))*step;
% Tdec=ceil(Dec/(J2*step))*step;
% J1=Acc/Tacc;
% J2=Dec/Tdec;
% 
% pa=1/Acc+1/Dec;
% pb=Tacc+Tdec;
% pc=-2*S;
% Vmax=roots([pa pb pc]);
% Vmax=Vmax(find(Vmax>0));
% 
% T1=ceil((Vmax/Acc-Tacc)/step)*step;
% T2=ceil((Vmax/Dec-Tdec)/step)*step;
% C1=Va*Tacc+Va*T1/2+Va*Tdec+Va*T2/2;
% Vmax=S/(Tacc+T1/2+Tdec+T2/2);
% %求解C2
% T2o=(Acc*Tacc-Dec*Tdec)/Dec;
% C2=Acc*Tacc^2+Acc*Tacc*Tdec+Acc*Tacc*T2o/2;
% %求解C2完成
% T3=0;
% T=2*(Tacc+Tdec)+T1+T2+T3;
% 
% t=0:step:T;
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% %v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% %d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%motioncase 3
% clear;%在C2和C3两个临界之间的情况：C2表示加速度恰好达到最大而减速度达到最大，C3表示加速度未达到最大而减速度恰好达到最大
% figure(3);
% J1=3;J2=2;
% step=0.001;
% %Va=6;   %机器允许的最大速度
% S=6;
% Acc=3;
% Dec=2;  %条件结束
% Tdec=ceil(Dec/(J2*step))*step;
% J2=Dec/Tdec;
% Tacc=real(roots([J1^2/(2*J2*Tdec) J1 J1*Tdec/2 0 -S]));
% Tacc=Tacc(4);
% Vmax=J1*Tacc^2;
% T1=0;T3=0;
% T2=ceil((J1*Tacc^2/J2*Tdec-Tdec)/step)*step;
% T=2*(Tacc+Tdec)+T1+T2+T3;
% %求解C3
% C3=J1*(Dec*Tdec/J1)^1.5+Dec*Tdec^2;
% %求解C3结束
% t=0:step:T;
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% %v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% %d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   C3=3.6330  有问题
% motioncase 4
% clear;
% figure(4);
% J1=3;J2=2;
% step=0.01;
% Va=3;   %机器允许的最大速度
% S=3;
% Tdec=ceil(((S/(J2*(1+(J2/J1)^0.5)))^(1/3))/step)*step;
% Tacc=ceil((J2/J1)^0.5*Tdec/step)*step;
% Acc=Tacc*J1;
% Dec=Tdec*J2;
% T1=0;T2=0;T3=0;
% T=2*(Tacc+Tdec)+T1+T2+T3;
% t=0:step:T;
% a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% motioncase 5  没问题
% clear;
% figure(5);
% J1=3;J2=2;
% step=0.01;
% V=3;   %机器允许的最大速度
% S=10;
% Acc=3;
% Dec=2;
% Tacc=ceil((V/J1)^0.5/step)*step;
% J1=V/Tacc;
% 
% T1=0;
% Tdec=ceil(Dec/J2/step)*step;
% T2=ceil((V/(Tdec*J2)-Tdec)/step)*step;
% J2=V/(Tdec*(T2+Tdec));
% J2=Tacc^2*J1/(Tdec*(T2+Tdec));
% C4=V*(Tacc+Tdec)+V*T2/3;
% T3=ceil((S/(Tacc^2*J1)-Tacc-Tdec-T2/2)/step)*step;
% J1=S/(Tacc^2*(T3+Tacc+Tdec+T2/2));
% J2=Tacc^2*J1/(Tdec*(T2+Tdec));
% T=2*(Tacc+Tdec)+T1+T2+T3;
% t=0:step:T;
% a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %motioncase 6  没问题
% clear;
% figure(6);
% J1=3;J2=2;
% step=0.01;
% Va=3;   %机器允许的最大速度
% S=7;
% Tacc=(Va/J1)^0.5;
% Tdec=roots([0.5*J1*Tacc^2  J1*Tacc^3-S  0.5*J1^2/J2*Tacc^4]);
% Tdec=Tdec(2);
% Acc=Tacc*J1;
% Dec=Tdec*J2;
% Tdec=ceil((Dec/J2)/step)*step;
% J2=Dec/Tdec;
% Tacc=real(roots([0.5*J1^2/J2/Tdec  J1  0.5*J1*Tdec 0 -S]));
% Tacc=Tacc(4);
% Acc=Tacc*J1;
% Dec=Tdec*J2;
% T1=0;T3=0;
% T2=J1*Tacc^2/(J2*Tdec)-Tdec;
% T=2*(Tacc+Tdec)+T1+T2+T3;
% t=0:step:T;
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% %v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% %d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%motioncase 7 
% clear;
% figure(7);
% J1=3;J2=2;
% step=0.01;
% S=2;
% Tdec=ceil(((S/(J2*(1+(J2/J1)^0.5)))^(1/3))/step)*step;
% Tacc=ceil(((J2/J1)^0.5*Tdec)/step)*step;
% Acc=Tacc*J1;
% Dec=Tdec*J2;
% T1=0;T2=0;T3=0;
% T=2*(Tacc+Tdec)+T1+T2+T3;
% C5=J1*(Dec*Tdec/J1)^1.5+Dec*Tdec^2;
% t=0:step:T;
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% %v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% %d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');

% Acc * Tacc > Dec * Tdec > V 
% % S > C6
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %case 8
% clear;
% figure(8);
% J1=3;J2=2;
% step=0.01;
% S=2;
% Acc=3;
% Dec=2;
% vmax=1;%允许的最大
% Tacc=ceil(((vmax/J1)^0.5)/step)*step;
% Tdec=ceil((vmax/J2)^0.5/step)*step;
%  Acc=Tacc*J1;
%  Dec=Tdec*J2;
% C6=vmax^1.5/((J1^-0.5)+(J2^-0.5));
% T3=(S-C6)/vmax;
% T=2*(Tacc+Tdec)+T3;
% T1=0;
% T2=0;
% Acc*Tacc
% Dec*Tdec
% t=0:step:T;
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
%    a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
%     a(i)=Acc-J1*(t(i)-Tacc-T1);
% elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
%     a(i)=0;
% elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T1-T3);
% elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
%     a(i)=-Dec;
% elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
%     a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
% end
% end
% subplot(311);
% plot(0:step:T-step,a);title('加速度');
% %v=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+a(i)*step;
%     end
%     v(k)=ss;
% end
% subplot(312);
% plot(0:step:T-step,v);title('速度');
% %d=[];
% for k=1:T/step
%     ss=0;
%     for i=1:k
%        ss=ss+v(i)*step;
%     end
%     d(k)=ss;
% end
% subplot(313);
% plot(0:step:T-step,d);title('位移');
%%%%%%%%%%%%%%
% case 9
% Acc * Tacc > Dec * Tdec > V 
% C6 > S
clear;
figure(9);
J1=3;J2=2;
step=0.01;
S=1;
Tdec=ceil(((S/(J2*(1+(J2/J1)^0.5)))^(1/3))/step)*step;
Tacc=ceil(((J2/J1)^0.5*Tdec)/step)*step;
S1=J1*Tacc^3+J2*Tdec^3;
J1=S1/(Tacc^3+Tacc^2*Tdec);
J2=J1*Tacc^2/Tdec^2;
T=2*(Tacc+Tdec);
t=0:step:T;
a=[];
for i=1:T/step
if(t(i)<=Tacc)
    a(i)=t(i)*J1;
elseif(Tacc<t(i)&&t(i)<=2*Tacc)
    a(i)=Tacc*J1-(t(i)-Tacc)*J1;
elseif(2*Tacc<t(i)&&t(i)<2*Tacc+Tdec)
    a(i)=-J2*(t(i)-2*Tacc);
elseif(2*Tacc+Tdec<=t(i))
    a(i)=-Tdec*J2+J2*(t(i)-2*Tacc-Tdec);
end
end
subplot(311);
plot(0:step:T-step,a);title('加速度');
v=[];
for k=1:T/step
    ss=0;
    for i=1:k
       ss=ss+a(i)*step;
    end
    v(k)=ss;
end
subplot(312);
plot(0:step:T-step,v);title('速度');
d=[];
for k=1:T/step
    ss=0;
    for i=1:k
       ss=ss+v(i)*step;
    end
    d(k)=ss;
end
subplot(313);
plot(0:step:T-step,d);title('位移');

