F= 8;
TM = @(t,x) L96(t,x,F);

ndim = 40;
nobs = 1000;
deltaobs = 0.05;

xrand = rand(ndim,1);
ttrans = linspace(0,100,1001);
% //TOCHANGE which to use, TF or TM?
[~,xtrans] = ode45(TM,ttrans,xrand);

xtrans = xtrans';
xattr = xtrans(:,end);
tend = nobs * deltaobs;
tobs = linspace(0,tend,nobs+1);
% 100 times finer for plotting.
trajfiner = 100;
ttraj = linspace(0,tend,nobs*trajfiner+1);

[~,truetraj] = ode45(TM,ttraj,xattr);
truestate = truetraj(1:trajfiner:end,:);

dlmwrite(strcat('Data/State','AllTraj.txt'),truestate)