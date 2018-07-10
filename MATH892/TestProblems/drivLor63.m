x0_ = [0,1,0];
x0 = 10*rand(1,3)-[10,10,10];
%x0 = [1,0,0];

b = 28; 
a = 10; 
c = 8/3;

xe1 = [sqrt(c*(b-1)),sqrt(c*(b-1)),b-1];
xe2 = [-sqrt(c*(b-1)),-sqrt(c*(b-1)),b-1];
xe3 = [0,0,0];

[t,x]=ode45(@FLor63,[0 20],x0);
figure 
hold on
plot3(x(:,1),x(:,2),x(:,3));
plot3(xe1(1),xe1(2),xe1(3),'b-*')
plot3(xe2(1),xe2(2),xe2(3),'r-*')
plot3(xe3(1),xe3(2),xe3(3),'k-*')
hold off
view(3)