data=dlmread('C:\Users\hqi\Desktop\logged_data.txt','',5);
subplot(411);
plot(data(:,1));title('wave1');
subplot(412);
plot(data(:,2));title('wave2');
subplot(413);
plot(data(:,3));title('wave3');
subplot(414);
plot(data(:,4));title('wave4');