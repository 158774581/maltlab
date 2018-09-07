% %结合在一起的总程序
% %motioncase 1
% clear;%速度达到最大，有匀速运动
% figure(1);
% J1=3;J2=2;
% step=0.01;
% Va=6;   %机器允许的最大速度
% S=22;
% Acc=3;
% Dec=2;
% Tdec=Dec/J2;
% Tacc=Acc/J1;
% T1=(Va-Acc*Tacc)/Acc;
% T2=(Va-Dec*Tdec)/Dec;
% C1=Va*Tacc+Va*T1/2+Va*Tdec+Va*T2/2;
% T3=(S-C1)/Va;
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
% %vmax=S/(Tacc+Tdec+0.5*T1+0.5*T2+T3)     //可以反推出vmax
%motioncase 2
% clear;%速度没有达到最大，无匀速运动。加（减）速度达到最大
% figure(2);
% J1=3;J2=2;
% step=0.01;
% Va=6;   %机器允许的最大速度
% S=18;
% Acc=3;
% Dec=2;
% Tdec=Dec/J2;
% Tacc=Acc/J1;
% pa=1/Acc+1/Dec;
% pb=Tacc+Tdec;
% pc=-2*S;
% Vmax=roots([pa pb pc]);
% Vmax=Vmax(find(Vmax>0));
% 
% T1=Vmax/Acc-Tacc;
% T2=Vmax/Dec-Tdec;
% C1=Va*Tacc+Va*T1/2+Va*Tdec+Va*T2/2;
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
% vmax=S/(Tacc+Tdec+0.5*T1+0.5*T2+T3)
% J11 = Vmax / (Tacc^2 + Tacc * T1)
%motioncase 3
clear;%在C2和C3两个临界之间的情况：C2表示加速度恰好达到最大而减速度达到最大，C3表示加速度未达到最大而减速度恰好达到最大
figure(3);
J1=3;J2=2;
step=0.01;
%Va=6;   %机器允许的最大速度
%S=6;
Acc=3;
Dec=2;  %条件结束
Tdec=Dec/J2;
Tacc=Acc/J1;
Vmax=J1*Tacc^2;
T1=0;T3=0;
T2=J1*Tacc^2/J2*Tdec-Tdec;
T=2*(Tacc+Tdec)+T1+T2+T3;
%求解C3
C3=J1*(Dec*Tdec/J1)^1.5+Dec*Tdec^2;
%求解C3结束
t=0:step:T;
%a=[];
for i=1:T/step
if(t(i)<=Tacc)
    a(i)=t(i)*J1;
elseif(Tacc<t(i)&&t(i)<=Tacc+T1)
   a(i)=Acc; %a(i)=Tacc*J1-(t(i)-Tacc)*J1;
elseif(Tacc+T1<t(i)&&t(i)<=2*Tacc+T1)
    a(i)=Acc-J1*(t(i)-Tacc-T1);
elseif(2*Tacc+T1<t(i)&&t(i)<=2*Tacc+T1+T3)
    a(i)=0;
elseif(2*Tacc+T1+T3<t(i)&&t(i)<=2*Tacc+T1+T3+Tdec)
    a(i)=-J2*(t(i)-2*Tacc-T1-T3);
elseif(2*Tacc+T1+T3+Tdec<t(i)&&t(i)<=2*Tacc+T1+T2+Tdec+T3)
    a(i)=-Dec;
elseif(2*Tacc+T1+T2+Tdec+T3<t(i))
    a(i)=-Dec+J2*(t(i)-(2*Tacc+T1+T2+Tdec+T3));   
end
end
subplot(311);
plot(0:step:T-step,a);title('加速度');
%v=[];
for k=1:T/step
    ss=0;
    for i=1:k
       ss=ss+a(i)*step;
    end
    v(k)=ss;
end
subplot(312);
plot(0:step:T-step,v);title('速度');
%d=[];
for k=1:T/step
    ss=0;
    for i=1:k
       ss=ss+v(i)*step;
    end
    d(k)=ss;
end
subplot(313);
plot(0:step:T-step,d);title('位移');
vmax=/(Tacc+Tdec+0.5*T1+0.5*T2+T3)
J11 = Vmax / (Tacc^2 + Tacc * T1)
%motioncase 4
% clear;
% figure(4);
% J1=3;J2=2;
% step=0.01;
% %Va=3;   %机器允许的最大速度
% S=10;
% Tdec=(S/(J2*(1+(J2/J1)^0.5)))^1/3;
% Tacc=(J2/J1)^0.5*Tdec;
% Acc=Tacc*J1;
% Dec=Tdec*J2;
% T1=0;T2=0;T3=0;
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
 %motioncase 5
% clear;
% figure(5);
% J1=3;J2=2;
% step=0.01;
% V=3;   %机器允许的最大速度
% S=10;
% Acc=3;
% Dec=2;
% Tacc=(V/J1)^0.5;
% T1=0;
% Tdec=Dec/J2;
% T2=(V-Tdec*Dec)/Dec;
% C4=V*(Tacc+Tdec)+V*T2/3;
% T3=(S-C4)/V;
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
% %motioncase 6
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
% clear;
% figure(7);
% J1=3;J2=2;
% step=0.01;
% S=6;
% Tdec=(S/(J2*(1+(J2/J1)^0.5)))^1/3;
% Tacc=(J2/J1)^0.5*Tdec;
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
% %case 8
% % Acc * Tacc > Dec * Tdec > V 
% % S > C6
% clear;
% figure(8);
% J1=3;J2=2;
% step=0.01;
% S=2.5;
% Acc=3;
% Dec=2;
% vmax=Dec^2/J2;%运动时的最大速度，不是允许的最大
% Tacc=(vmax/J1)^0.5;
% Tdec=(vmax/J2)^0.5;
% C6=vmax^1.5/(1/J1^0.5+1/J2^0.5);
% T3=(S-C6)/vmax;
% T=2*(Tacc+Tdec)+T3;
% t=0:step:T;
% 
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=2*Tacc)
%     a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(2*Tacc<t(i)&&t(i)<2*Tacc+T3)
%     a(i)=0;
% elseif(2*Tacc+T3<=t(i)&&t(i)<2*Tacc+T3+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc-T3);
% elseif(2*Tacc+T3+Tdec<=t(i))
%     a(i)=-Tdec*J2+J2.*(t(i)-2*Tacc-Tdec-T3);
% end
% end
%  subplot(311);
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
% %case 9
% % Acc * Tacc > Dec * Tdec > V 
% % C6 > S
% clear;
% figure(9);
% J1=3;J2=2;
% step=0.01;
% S=3;
% Tdec=(S/(J2*(1+(J2/J1)^0.5)))^1/3;
% Tacc=(J2/J1)^0.5*Tdec;
% T=2*(Tacc+Tdec);
% t=0:step:T;
% %a=[];
% for i=1:T/step
% if(t(i)<=Tacc)
%     a(i)=t(i)*J1;
% elseif(Tacc<t(i)&&t(i)<=2*Tacc)
%     a(i)=Tacc*J1-(t(i)-Tacc)*J1;
% elseif(2*Tacc<t(i)&&t(i)<2*Tacc+Tdec)
%     a(i)=-J2*(t(i)-2*Tacc);
% elseif(2*Tacc+Tdec<=t(i))
%     a(i)=-Tdec*J2+J2*(t(i)-2*Tacc-Tdec);
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
% 
