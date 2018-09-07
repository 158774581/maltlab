 syms t m2 l1 l2 et1 et2 t1 t2 y;
[et1,et2]=dsolve('t2=m2*l1*l2*(cos(et2)*D2et1)','t1=m2*l2^2*(D2et1+D2et2)','et1(0)=0','Det1(0)=1','et2(0)=0','Det2(0)=1')



