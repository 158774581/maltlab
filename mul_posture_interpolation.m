R=zeros(3,3,10);
a=zeros(3,3,10);
b=zeros(3,3,10);
c=zeros(3,3,10);
s=zeros(3,3);
vec_s=zeros(3,1);
vec_u=zeros(3,1);
vec_t=zeros(3,1);
t=zeros(3,3);
u=zeros(3,3);
R_tmp=zeros(3,3,10);
r=zeros(3,3,10);
r_j=zeros(3,3,10);
alpha=zeros(10,1);
fai=zeros(10,1);
alpha(:)=pi/3/10:pi/3/10:pi/3;
belta=zeros(10,1);
belta(:)=pi/2/10:pi/2/10:pi/2;
gamma=zeros(10,1);
gamma(:)=pi/10:pi/10:pi;
%calculate R1...R10
for i=1:10
    R(:,:,i)=[cos(alpha(i))*cos(belta(i)) cos(alpha(i))*sin(belta(i))*sin(gamma(i))-sin(alpha(i))*cos(gamma(i)) cos(alpha(i))*sin(belta(i))*cos(gamma(i))+sin(alpha(i))*sin(gamma(i))
              sin(alpha(i))*cos(belta(i)) sin(alpha(i))*sin(belta(i))*sin(gamma(i))+cos(alpha(i))*cos(gamma(i)) cos(alpha(i))*sin(belta(i))*cos(gamma(i))-cos(alpha(i))*sin(gamma(i))
              -sin(belta(i)) cos(belta(i))*sin(gamma(i)) cos(belta(i))*cos(gamma(i))];
end
for i=2:10
    R_tmp(:,:,i)=R(:,:,i-1)'*R(:,:,i);
     fai(i)=acos((R_tmp(1,1,i)+R_tmp(2,2,i)+R_tmp(3,3,i)-1)/2);
          r_j(:,:,i)=( R_tmp(:,:,i)- R_tmp(:,:,i)')/2/sin(fai(i));
          r(:,:,i)=r_j(:,:,i)*fai(i);
end
w1=1;w2=0.5;
c(:,:,2)=[0 0 0; 0 0 -w1; 0 w1 0];
b(:,:,2)=[0 w2 0; w2 0 0; 0 0 0];
a(:,:,2)=r(:,:,2)-b(:,:,2)-c(:,:,2);
for i=3:10
    s(:,:)=r(:,:,i);
    vec_s=[s(3,2);s(1,3);s(2,1)];
    norm=(s(2,3)^2+s(3,1)^2+s(2,1)^2)^0.5;
    t=3*a(:,:,i-1)+2*b(:,:,i-1)+c(:,:,i-1);
    vec_t=[t(3,2);t(1,3);t(2,1)];
    u=6*a(:,:,i-1)+2*b(:,:,i-1);
    vec_u=[u(3,2);u(1,3);u(2,1)];
    temp1=zeros(3,3);
    temp2=zeros(3,3);
    temp3=zeros(3,3);
    temp0=zeros(3,1);
    temp0=(1-cos(norm))/norm^2*cross(vec_s,vec_u);
    temp1=[0 -temp0(3) temp0(2); temp0(3) 0 -temp0(1);temp0(2) temp0(1) 0];
    temp0=(3*sin(norm)-norm*cos(norm)-2*norm)*cross(vec_s,vec_s-vec_t)/norm^5;
    temp2=[0 -temp0(3) temp0(2); temp0(3) 0 -temp0(1);temp0(2) temp0(1) 0]*s'*t;
    temp0=(norm-sin(norm))/norm^3*(cross(vec_t,vec_s-vec_t)+cross(vec_s,cross(vec_s,vec_u)));
    temp3=[0 -temp0(3) temp0(2); temp0(3) 0 -temp0(1);temp0(2) temp0(1) 0];
    c(:,:,i)=([1 0 0;0 1 0;0 0 1]-(1-cos(norm))/norm^2*s+(norm-sin(norm))/norm^3*s*s)*t;
    b(:,:,i)=0.5*(u-s'*t/norm^5*(2*cos(norm)+norm*sin(norm)-2)*(s-t)-temp1+temp2+temp3);
    a(:,:,i)=s- b(:,:,i)-c(:,:,i);
end
for i=2:10
%r=a+b+c
temp1=a(:,:,i)+b(:,:,i)+c(:,:,i);
norm=(temp1(2,3)^2+temp1(3,1)^2+temp1(2,1)^2)^0.5;
%Ri
temp2=[1 0 0;0 1 0;0 0 1]+sin(norm)/norm*temp1+(1-cos(norm))/norm^2*temp1*temp1;
%R(t)=Ri-1*R
i
temp3= R(:,:,i-1)*temp2;
belta1=atan2(-temp3(3,1),(temp3(1,1)^2+temp3(2,1)^2)^0.5);
belta2=atan2(-R(3,1,i),(R(1,1,i)^2+R(2,1,i)^2)^0.5);
end
% 2%_3%µÄÎó²î¡£

          