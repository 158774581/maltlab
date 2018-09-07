a=[0 1 2];    %the length of the connecting rod
alpha=[0 0 0];    
d=[0 0 0];
theata=[10*pi/180 20*pi/180 30*pi/180];
T=zeros(4,4,3);    %i corrresponding to 01T
for i=1:3
    T(:,:,i)=[cos(theata(i)) -sin(theata(i)) 0 a(i)
        sin(theata(i))*cos(alpha(i)) cos(theata(i))*cos(alpha(i)) -sin(alpha(i)) -sin(alpha(i))*d(i)
         sin(theata(i))*sin(alpha(i)) cos(theata(i))*sin(alpha(i))  cos(alpha(i)) cos(alpha(i))*d(i)
         0 0 0 1]
end
T0=T(:,:,1)*T(:,:,2)*T(:,:,3)    %positive solution
x=T0(1,4);y=T0(2,4);
l1=1;l2=2;
theta2=acos((x^2+y^2-l1^2-l2^2)/2/l1/l2)    %reverse solution of theata2

