function DataGeneration(rngs)
%% Parameters allowed for change

F = 8;

nobs = 200;
deltaobs = 0.05;

ndim = 40;

dobs = 20;

nens = 200;

%% Other parameters

rng(rngs);

% True model
TM = @(t,x) L96(t,x,F);

% Forecast models
TF = @(t,x) L96(t,x,F);

% Note: They could be different.
% That was something I was working on.
% But for now we just take the same.

% Observation operators
% H = [eye(dobs),zeros(dobs,ndim-dobs)];
h = @(x) obs_40_20(x);

% stupid inverse case
% h = @(x) 

% Seems like nonlinear not working, try linear
% h = @(x) x(1:20,:);
% H = [eye(dobs),zeros(dobs,ndim-dobs)];

% Magic number that works very well
sigmaobs = 0.00;

%% 

xrand = rand(ndim,1);
ttrans = linspace(0,100,1001);
[~,xtrans] = ode45(TM,ttrans,xrand);

xtrans = xtrans';
xattr = xtrans(:,end);


tend = nobs * deltaobs;
tobs = linspace(0,tend,nobs+1);
% 100 times finer for plotting.
trajfiner = 100;
ttraj = linspace(0,tend,nobs*trajfiner+1);

[~,truetraj] = ode45(TM,ttraj,xattr);
truetraj = truetraj';
truestate = truetraj(:,1:trajfiner:end);

%% Generate obervations

% Assume diagonal covariance
% Change from linear to full
% trueobs = H*truestate;
trueobs = h(truestate);
% robs = sigmaobs^2;
% Robsmat = eye(dobs)*robs;
% yobs = trueobs+mvnrnd(zeros(dobs,1),Robsmat,nobs+1)';

%%

dlmwrite(sprintf('Data/obs%d.txt',rngs),trueobs)
dlmwrite(sprintf('Data/states%d.txt',rngs),truestate)

end