%% FTHS
dim=9;
x0 = randn(dim,1);

[t,x]=ode45(@(t,u)FTHS(t,u)',[0 20],x0);
plot3(x(:,1),x(:,2),x(:,3));

%% KSE
dim = 3;
x0 = randn(1,dim);
x0 = zeros(1,dim);
x0(1)=1;

[t,x]=ode113(@(t,u)kse(t,u,32)',[0 2e-1],x0);
plot3(x(:,1),x(:,2),x(:,3));