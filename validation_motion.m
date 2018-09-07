spd_max=6;
acc_max=3;
dec_max=1;
Tacc=1;
Tdec=1;
Jacc=3;
Jdec=1;
pos=8;
c0=(2 * Jdec*Tdec) / Jacc;
c1 = (Jdec*Tdec*Tdec) /Jacc;
c2 = 0;
c3 = (-2 * pos*Jdec*Tdec) / (Jacc*Jacc);
tmp  = 0.25 *c0;
tmp1 =c0 * c0;
p = c1 - 0.375* tmp1;
q = (0.125*tmp1-0.5*c1)* c0 + c2;
r = (-0.01171875* tmp1*tmp1) + c3 - (0.25*c0*c2) + (0.0625*tmp1*c1);
    c0 = p;
	c1 = 0.25*p*p-r;
	c2 = -0.125*q*q;%三次方程系数3个;
    
    %化为y^3+p*y+q
    tmp = 0.3333333333333333 * c0 ;
    p = c1 - tmp*c0;
	q = 0.0740740740740741*c0 *c0*c0 - tmp*c1 + c2;
    discriminant = 0.037037037037037 * p * p * p + 0.25 * q * q;
    u3 = -0.5*q + discriminant;
		v3 = -0.5*q - discriminant;
        u = u3^0.33333333333333333;
        v = -(-v3)^0.33333333333333333;
    t1 = u+v;
    r0 = t1-tmp;%CUBIC_ONE_REAL_TWO_COMPLEX_ROOTS  得到3.2.1的根
    m=r0;
    c0 = (2*m)^0.5;
	c1 = 0.5*p+m-(0.5*q)/(c0);
  bot=  roots([1 c0 c1])
discriminant1 = 0.25 * c0 *c0 - c1;
c0 = -(2*m)^0.5;
c1 = 0.5*p+m+(0.5*q)/(c0);
    discriminant2 = 0.25 * c0 * c0 - c1 ;
    
    
    
    
    
    
    
