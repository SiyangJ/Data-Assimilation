x0 = randn(1,40);
x0 = zeros(1,40);
x0(1)=1;

[t,x]=ode45(@FLor95,[0 20],x0);
plot3(x(:,1),x(:,2),x(:,3));
