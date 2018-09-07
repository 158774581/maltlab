dtheata=zeros(1000,3); 
theata=zeros(1000,3);
X=[0.2 -0.3 -0.2]';    %dx/dt dy/dt w
theata(1,:)=[10*pi/180 20*pi/180 30*pi/180];
L1=4;L2=3;L3=2;
jacob=zeros(1000,3,3);
jacobi=zeros(3,3);
x=zeros(1000); 
y=zeros(1000);
J=zeros(1000);  %det(jacobi)
T=zeros(1000,3);    % joint torque
for i=1:50
jacob(i,:,:)=[-L1*sin(theata(i,1))-L2*sin(theata(i,1)+theata(i,2))-L3*sin(theata(i,1)+theata(i,2)+theata(i,3)) ,-L2*sin(theata(i,1)+theata(i,2))-L3*sin(theata(i,1)+theata(i,2)+theata(i,3)),-L3*sin(theata(i,1)+theata(i,2)+theata(i,3))
L1*cos(theata(i,1))+L2*cos(theata(i,1)+theata(i,2))+L3*cos(theata(i,1)+theata(i,2)+theata(i,3)) ,L2*cos(theata(i,1)+theata(i,2))+L3*cos(theata(i,1)+theata(i,2)+theata(i,3)),L3*cos(theata(i,1)+theata(i,2)+theata(i,3))
1 1 1];
x(i)=L1*cos(theata(i,1))+L2*cos(theata(i,1)+theata(i,2))+L3*cos(theata(i,1)+theata(i,2)+theata(i,3));
y(i)=L1*sin(theata(i,1))+L2*sin(theata(i,1)+theata(i,2))+L3*sin(theata(i,1)+theata(i,2)+theata(i,3)) ;
jacobi(:,:)=jacob(i,:,:);
T(i,:)=jacobi'*[1 2 3]';
J(i)=det(jacobi);
dtheata(i,:)=((jacobi)\X);
theata(i+1,:)=theata(i,:)+dtheata(i,:)*0.1;
end
figure(1);
subplot(411);
plot(0:0.1:5,dtheata(1:51,1));
subplot(412);
plot(0:0.1:5,theata(1:51,1));
subplot(413);
plot(0:0.1:5,x(1:51));
subplot(414);
plot(0:0.1:5,y(1:51));
figure(2);
subplot(411);
plot(0:0.1:5,J(1:51));
subplot(412)
plot(0:0.1:5,T(1:51,3));

