%% Parameters allowed for change

F = 8;
A = 0*F;
B = 0*F;

nobs = 50;
deltaobs = 0.1;

ndim = 40;

dobs = 10;

nens = 200;

inflmu = 1e-4;

%% Other parameters

rng(1);


sinus = sin(2*pi/ndim*(0:(ndim-1))');
zeta = A*sinus;
xi = B*sinus;

TM = @(t,x) L96ME(t,x,F,xi,zeta);
TF = @(t,x) L96(t,x,F);

H = [eye(dobs),zeros(dobs,ndim-dobs)];

%% Find an initial condition on the attractor

xrand = rand(ndim,1);
ttrans = linspace(0,100,1001);
% //TOCHANGE which to use, TF or TM?
[~,xtrans] = ode45(TM,ttrans,xrand);

xtrans = xtrans';
xattr = xtrans(:,end);

%% Plot a long trajectory

tlong = linspace(0,100,1e4);
% //TOCHANGE 
[~,xlong] = ode45(TM,tlong,xattr);
xlong = xlong';

%sigmaobs = (max(xlong(1,:))-min(xlong(1,:)))/50;
%sigmainit = sigmaobs*10;
sigmaobs = 0.09;

%% Generate truth 

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
trueobs = H*truestate;
robs = sigmaobs^2;
Robsmat = eye(dobs)*robs;
yobs = trueobs+mvnrnd(zeros(dobs,1),Robsmat,nobs+1)';








%% Generate ensemble

% Recommended choice by Baek 2006
sigmainit = 1.3;

xfm = zeros(ndim,nobs+1);
xam = zeros(ndim,nobs+1);
xfcov = zeros(ndim,ndim,nobs+1);
xacov = zeros(ndim,ndim,nobs+1);

xmf = zeros(ndim,nobs+1);
xcovf = zeros(ndim,ndim,nobs+1);

% Initial mean and covariance
xam(:,1) = xattr + sigmainit*randn(ndim,1);
xacov(:,:,1) = sigmainit^2*eye(ndim);

xmf(:,1) = xam(:,1);
xcovf(:,:,1) = xacov(:,:,1);

% Freerun ensemble
xensf = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (perfect model assumption)
xens = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

%% Perform EnKF

for i=2:nobs+1
    
    % Forecast
    % Evolve the ensemble
    for j = 1:nens   
        [~,xf] = ode45(TF,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xens(:,j) = xf; 
        
        [~,xf] = ode45(TF,tobs(i-1:i),xensf(:,j));
        xf = xf(end,:)';
        xensf(:,j) = xf;
    end
    
    % Calculate the stats of ensemble for xf
    xfm(:,i) = mean(xens,2);
    xfcov(:,:,i) = cov(xens');
    
    xmf(:,i) = mean(xensf,2);
    xcovf(:,:,i) = cov(xens');
    
    % Update
    % Kalman gain
    
    % Enhanced variance inflation
    % Ott et al. 2004
    pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
    pfht = pf*H';
    
    kgain = pfht*(H*pfht+Robsmat)^-1;
    
    for j=1:nens
        
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xens(:,j) = xens(:,j)+kgain*(pertob-H*xens(:,j));
        
    end
    
    % Calculate the stats of ensemble for xa
    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
end


%% Other stuff

figure
for i=1:4
    subplot(1,4,i)
    hold on
    plot(ttraj,truetraj(i,:),'b-','LineWidth',4)
    plot(tobs,xam(i,:),'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    plot(tobs,xfm(i,:),'r-s','LineWidth',2,'MarkerSize',8)
    plot(tobs,xmf(i,:),'g--*')
    hold off
end