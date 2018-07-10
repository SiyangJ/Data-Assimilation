function BM = EnKFfun_ME2_v2(mu,gamma,fname)
%ENKFFUN_V1 Summary of this function goes here
%   Detailed explanation goes here

% Used for Model Error 2
% Save the local workspace for inspection
% Doesn't specifically return the benchmarks

% All stuff about Bias Model 2 are commented out
% Baek 2006 suggests Bias Model 3 might beat BM 1
% So we keep BM 3

DEBUG__ = true;

%% Parameters allowed for change

fprintf('Start: mu=%f\n',...
    mu);

% For Model Error 2,
% i.e. state-dependent model error 
% -gamma * x_i^2
% Baek 2006 suggests use F = 10 to maintain chaotic behavior
% For zeta and xi, we just use 0
F = 10;
%gamma = 0.05

% For Type A error, 100 observations are enough for error to settle
% Model Error 2 is similar, so we might as well use 100

% For testing purposes, we use 800

nobs = 500;
deltaobs = 0.05;

ndim = 40;

dobs = 40;

% Seems a reasonable choice
% Could try others like 60
nens = 40;

inflmu = mu;

%% Other parameters

rng(1);

% True model
TM = @(t,x) L96ME2(t,x,F,0,0,gamma);
% Forecast models
TF = @(t,x) L96(t,x,F);
TB = @(t,x,ga) L96ME2(t,x,F,0,0,ga);

% Observation operators
H = [eye(dobs),zeros(dobs,ndim-dobs)];

Hb1 = [eye(dobs),zeros(dobs,ndim+1-dobs)];

% For bias model 3, [x;zeta;xi]

sigmaobs = 0.09;

% Recommended choice by Baek 2006
sigmainit = 1.3;

% Need to test
sigmagamma = 1.0;

%% Find an initial condition on the attractor

xrand = rand(ndim,1);
ttrans = linspace(0,100,1001);
[~,xtrans] = ode45(TM,ttrans,xrand);

xtrans = xtrans';
xattr = xtrans(:,end);

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

% No bias handling
xfm = zeros(ndim,nobs+1);
xam = zeros(ndim,nobs+1);
xfcov = zeros(ndim,ndim,nobs+1);
xacov = zeros(ndim,ndim,nobs+1);

% Free run
xmf = zeros(ndim,nobs+1);
xcovf = zeros(ndim,ndim,nobs+1);

% Bias model I
xfmb1 = zeros(ndim+1,nobs+1);
xamb1 = zeros(ndim+1,nobs+1);
xfcovb1 = zeros(ndim+1,ndim+1,nobs+1);
xacovb1 = zeros(ndim+1,ndim+1,nobs+1);

% Initial mean and covariance
% Need to test initial mean
xam(:,1) = xattr + sigmainit*randn(ndim,1);
xacov(:,:,1) = sigmainit^2*eye(ndim);

xmf(:,1) = xam(:,1);
xcovf(:,:,1) = xacov(:,:,1);

% For bias models
% Initial estimate for zeta and xi are just zero
xamb1(1:ndim,1) = xam(:,1);
xacovb1(1:ndim,1:ndim,1) = xacov(:,:,1);
xacovb1(ndim+1,ndim+1,1) = sigmagamma^2;

% Freerun ensemble
xensf = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (perfect model assumption)
xens = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (bias model 1)
xensb1 = mvnrnd(xamb1(:,1),xacovb1(:,:,1),nens)';

%% Perform EnKF

for i=2:nobs+1
    
    % Forecast
    % Evolve the ensemble
    
    % Free run and Perfect Model
    for j = 1:nens          
        [~,xf] = ode45(TF,tobs(i-1:i),xensf(:,j));
        xf = xf(end,:)';
        xensf(:,j) = xf;
        
        [~,xf] = ode45(TF,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xens(:,j) = xf;  
        
        gtemp = xensb1(end,j);
        
        [~,xf] = ode45(@(t,x)TB(t,x,gtemp),...
            tobs(i-1:i),xensb1(1:ndim,j));
        xf = xf(end,:)';
        xensb1(1:ndim,j) = xf;
    end
    
    % Calculate the stats of ensemble for xf
   
    xmf(:,i) = mean(xensf,2);
    xcovf(:,:,i) = cov(xens');
    
    xfm(:,i) = mean(xens,2);
    xfcov(:,:,i) = cov(xens');
    
    xfmb1(:,i) = mean(xensb1,2);
    xfcovb1(:,:,i) = cov(xensb1');
    
    % Update
    
    % Kalman gain
    
    % Perfect Model   
    % Enhanced variance inflation
    % Ott et al. 2004
    pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
    pfht = pf*H';   
    kgain = pfht*(H*pfht+Robsmat)^-1;   
    
    % Bias Model 1   
    pf = xfcovb1(:,:,i)+inflmu*trace(xfcovb1(:,:,i))/(ndim+1)*eye(ndim+1);
    pfht = pf*Hb1';   
    kgainb1 = pfht*(Hb1*pfht+Robsmat)^-1;
    
    % Update Ensemble members
    
    for j=1:nens     
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xens(:,j) = xens(:,j)+kgain*(pertob-H*xens(:,j)); 
        xensb1(:,j) = xensb1(:,j)+kgainb1*(pertob-Hb1*xensb1(:,j));            
    end  
    
    % Calculate the stats of ensemble for xa    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
    
    xamb1(:,i) = mean(xensb1,2);
    xacovb1(:,:,i) = cov(xensb1');
    
    if DEBUG__
        fprintf('/')
        if mod(i,50)==0
            fprintf('\n')
        end
    end
    
end

%% Calculate benchmarks
sumrmsap = 0;
sumrmsab1 = 0;

tspan = floor(nobs/2)+1:nobs+1;
T = length(tspan);

for t = tspan
    
    sumrmsap = sumrmsap + norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsab1 = sumrmsab1 + norm(xamb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
end

eap = sumrmsap/T;
eab1 = sumrmsab1/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

sumrmsp = 0;
sumrmsb1 = 0;

for t = tspan
    
    sumrmsp = sumrmsp + norm(xfm(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb1 = sumrmsb1 + norm(xfmb1(1:ndim,t)-truestate(:,t))/sqrt(ndim); 
end

ep = sumrmsp/T;
eb1 = sumrmsb1/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gammab1 = mean(xfmb1(end,tspan),2);
gammab1cov = diag(cov(xfmb1(end,tspan)')).^0.5;

BM = {inflmu,gamma,deltaobs,nens,...
     eap,eab1,...
     ep,eb1,...
     gammab1,gammab1cov};

save(fname)

fprintf('End: mu=%f\n',...
    mu);

end

