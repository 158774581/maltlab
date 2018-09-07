syms t;syms b;syms y;
y=(1-t)^5+5*t*(1-t)^4+9*t^2*(1-t)^3;
y=expand(y)
dy=diff(y,t,1);
dy





