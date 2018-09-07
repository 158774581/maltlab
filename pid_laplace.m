syms s;syms t;
g1=4.375/(0.0006*s+1.4);
J=0.00844+1/200^2;
C=0.00013+0.5/200^2;
g2=1/(J*s+C);
oo=g1*g2/(1+0.00867*g1*g2);
square=5*(1+2*s+1/(3*s))*12*oo*(1/(200*s));
transferfunc=square/(1+square);
in=40*t;
IN=laplace(in)
transferfunc*IN
u=ilaplace(transferfunc*IN);
u=vpa(u,6)
ezplot(u);
% i=0;
% for t=0:0.01:2
%     i=i+1;
% u(i)=282.392*exp(-251.352*t) - 282.399*exp(-2081.5*t) + 0.0067323*exp(-0.246595*t)*cos(0.32223*t) + 0.00107128*exp(-0.246595*t)*sin(0.32223*t);
% % tt=0:0.1:2;
% % z=subs(u,t,tt);
% end
% plot(u);