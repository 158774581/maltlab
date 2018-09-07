tacc=1;
T1=2;
T2=2;
T3=2;
tdec=1;
J=2;

for i=1:n
if(t(i)<tacc||t(i)>(2*tacc+T1+T3+T2+tdec))
    jerk(i)=J;
else
    if((t(i)<(T1+2*tacc)&&t(i)>(T1+tacc))||t(i)>(T1+2*tacc+T3)&&t(i)<(T1+2*tacc+T3+tdec));
        jerk(i)=-J;
    else
        jerk=0;
    end      
end
end
plot(jerk,t);