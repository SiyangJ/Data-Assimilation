function BM = EnKFfun_ME2_v1(mu,gamma,fname)
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
A = 0*F;
B = 0*F;
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

% Different ensemble size for different models
nensarr = [2,2,3]*nens;

inflmu = mu;

%% Other parameters

rng(1);

sinus = sin(2*pi/ndim*(0:(ndim-1))');
zeta = A*sinus;
xi = B*sinus;

% True model
TM = @(t,x) L96ME2(t,x,F,xi,zeta,gamma);
% Forecast models
TF = @(t,x) L96(t,x,F);

% Observation operators
H = [eye(dobs),zeros(dobs,ndim-dobs)];

Hb1 = [eye(dobs),zeros(dobs,2*ndim-dobs)];
%Hb2 = [eye(dobs),zeros(dobs,ndim-dobs),eye(dobs),zeros(dobs,ndim-dobs)];
Hb3 = [eye(dobs),zeros(dobs,2*ndim-dobs),eye(dobs),zeros(dobs,ndim-dobs)];

% For bias model 3, [x;zeta;xi]

sigmaobs = 0.09;

% Recommended choice by Baek 2006
sigmainit = 1.3;

% Need to test
sigmazeta = 1.0;

% Need to test
sigmaxi = 1.0;

%% Find an initial condition on the attractor

xrand = rand(ndim,1);
ttrans = linspace(0,100,1001);
% //TOCHANGE which to use, TF or TM?
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
xfmb1 = zeros(2*ndim,nobs+1);
xamb1 = zeros(2*ndim,nobs+1);
xfcovb1 = zeros(2*ndim,2*ndim,nobs+1);
xacovb1 = zeros(2*ndim,2*ndim,nobs+1);

% Bias model II
% xfmb2 = zeros(2*ndim,nobs+1);
% xamb2 = zeros(2*ndim,nobs+1);
% xfcovb2 = zeros(2*ndim,2*ndim,nobs+1);
% xacovb2 = zeros(2*ndim,2*ndim,nobs+1);

% Bias model III
xfmb3 = zeros(3*ndim,nobs+1);
xamb3 = zeros(3*ndim,nobs+1);
xfcovb3 = zeros(3*ndim,3*ndim,nobs+1);
xacovb3 = zeros(3*ndim,3*ndim,nobs+1);

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
xacovb1(ndim+1:end,ndim+1:end,1) = sigmazeta^2*eye(ndim);

% xamb2(1:ndim,1) = xam(:,1);
% xacovb2(1:ndim,1:ndim,1) = xacov(:,:,1);
% xacovb2(ndim+1:end,ndim+1:end,1) = sigmaxi^2*eye(ndim);

xamb3(1:ndim,1) = xam(:,1);
xacovb3(1:ndim,1:ndim,1) = xacov(:,:,1);
xacovb3(ndim+1:2*ndim,ndim+1:2*ndim,1) = sigmazeta^2*eye(ndim);
xacovb3(2*ndim+1:end,2*ndim+1:end,1) = sigmaxi^2*eye(ndim);

% Freerun ensemble
xensf = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (perfect model assumption)
xens = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (bias model 1)
xensb1 = mvnrnd(xamb1(:,1),xacovb1(:,:,1),nensarr(1))';

% Initial ensemble (bias model 2)
% xensb2 = mvnrnd(xamb2(:,1),xacovb2(:,:,1),nensarr(2))';

% Initial ensemble (bias model 3)
xensb3 = mvnrnd(xamb3(:,1),xacovb3(:,:,1),nensarr(3))';

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
    end
    
    % Bias Model 1
    for j=1:nensarr(1)        
        [~,xf] = ode45(TF,tobs(i-1:i),xensb1(1:ndim,j));
        xf = xf(end,:)';
        xensb1(1:ndim,j) = xf + xensb1(ndim+1:end,j); 
    end
    
    % Bias Model 2
    % for j=1:nensarr(2)       
    %     [~,xf] = ode45(TF,tobs(i-1:i),xensb2(1:ndim,j));
    %     xf = xf(end,:)';
    %     xensb2(1:ndim,j) = xf;  
    % end
    
    % Bias Model 3
    for j=1:nensarr(3)   
        [~,xf] = ode45(TF,tobs(i-1:i),xensb3(1:ndim,j));
        xf = xf(end,:)';
        xensb3(1:ndim,j) = xf + xensb3(ndim+1:2*ndim,j);        
    end
    
    % Calculate the stats of ensemble for xf
   
    xmf(:,i) = mean(xensf,2);
    xcovf(:,:,i) = cov(xens');
    
    xfm(:,i) = mean(xens,2);
    xfcov(:,:,i) = cov(xens');
    
    xfmb1(:,i) = mean(xensb1,2);
    xfcovb1(:,:,i) = cov(xensb1');
    
    % xfmb2(:,i) = mean(xensb2,2);
    % xfcovb2(:,:,i) = cov(xensb2');
    
    xfmb3(:,i) = mean(xensb3,2);
    xfcovb3(:,:,i) = cov(xensb3');
    
    % Update
    
    % Kalman gain
    
    % Perfect Model   
    % Enhanced variance inflation
    % Ott et al. 2004
    pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
    pfht = pf*H';   
    kgain = pfht*(H*pfht+Robsmat)^-1;   
    
    % Bias Model 1   
    pf = xfcovb1(:,:,i)+inflmu*trace(xfcovb1(:,:,i))/2/ndim*eye(2*ndim);
    pfht = pf*Hb1';   
    kgainb1 = pfht*(Hb1*pfht+Robsmat)^-1;   
    
    % Bias Model 2
    % pf = xfcovb2(:,:,i)+inflmu*trace(xfcovb2(:,:,i))/2/ndim*eye(2*ndim);
    % pfht = pf*Hb2';   
    % kgainb2 = pfht*(Hb2*pfht+Robsmat)^-1;    
    
    % Bias Model 3
    pf = xfcovb3(:,:,i)+inflmu*trace(xfcovb3(:,:,i))/3/ndim*eye(3*ndim);
    pfht = pf*Hb3';   
    kgainb3 = pfht*(Hb3*pfht+Robsmat)^-1;   
    
    % Update Ensemble members
    
    for j=1:nens     
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xens(:,j) = xens(:,j)+kgain*(pertob-H*xens(:,j));        
    end  
    
    for j=1:nensarr(1)   
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xensb1(:,j) = xensb1(:,j)+kgainb1*(pertob-Hb1*xensb1(:,j));            
    end
    
    % for j=1:nensarr(2)    
    %     pertob = mvnrnd(yobs(:,i),Robsmat,1)';      
    %     xensb2(:,j) = xensb2(:,j)+kgainb2*(pertob-Hb2*xensb2(:,j));       
    % end
    
    for j=1:nensarr(3)   
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';     
        xensb3(:,j) = xensb3(:,j)+kgainb3*(pertob-Hb3*xensb3(:,j));       
    end
    
    % Calculate the stats of ensemble for xa    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
    
    xamb1(:,i) = mean(xensb1,2);
    xacovb1(:,:,i) = cov(xensb1');
    
    % xamb2(:,i) = mean(xensb2,2);
    % xacovb2(:,:,i) = cov(xensb2');
    
    xamb3(:,i) = mean(xensb3,2);
    xacovb3(:,:,i) = cov(xensb3');
    
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
% sumrmsab2 = 0;
sumrmsab3 = 0;

tspan = floor(nobs/2)+1:nobs+1;
T = length(tspan);

for t = tspan
    
    sumrmsap = sumrmsap + norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsab1 = sumrmsab1 + norm(xamb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    % sumrmsab2 = sumrmsab2 + norm(xamb2(1:ndim,t)+xamb2(ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    sumrmsab3 = sumrmsab3 + norm(xamb3(1:ndim,t)+xamb3(2*ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    
end

eap = sumrmsap/T;
eab1 = sumrmsab1/T;
% eab2 = sumrmsab2/T;
eab3 = sumrmsab3/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

sumrmsp = 0;
sumrmsb1 = 0;
% sumrmsb2 = 0;
sumrmsb3 = 0;

for t = tspan
    
    sumrmsp = sumrmsp + norm(xfm(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb1 = sumrmsb1 + norm(xfmb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    % sumrmsb2 = sumrmsb2 + norm(xfmb2(1:ndim,t)+xfmb2(ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb3 = sumrmsb3 + norm(xfmb3(1:ndim,t)+xfmb3(2*ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    
end

ep = sumrmsp/T;
eb1 = sumrmsb1/T;
% eb2 = sumrmsb2/T;
eb3 = sumrmsb3/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zetab1 = mean(xfmb1(ndim+1:end,tspan),2);
zetab1cov = diag(cov(xfmb1(ndim+1:end,tspan)')).^0.5;
% xib2 = mean(xfmb2(ndim+1:end,tspan),2);
% xib2cov = diag(cov(xfmb2(ndim+1:end,tspan)')).^0.5;
zetab3 = mean(xfmb3(ndim+1:2*ndim,tspan),2);
zetab3cov = diag(cov(xfmb3(ndim+1:2*ndim,tspan)')).^0.5;
xib3 = mean(xfmb3(2*ndim+1:end,tspan),2);
xib3cov = diag(cov(xfmb3(2*ndim+1:end,tspan)')).^0.5;

BM = {inflmu,A,B,gamma,deltaobs,nens,...
     eap,eab1,eab3,...
     ep,eb1,eb3,...
     zetab1,zetab1cov,zetab3,zetab3cov,xib3,xib3cov};

save(fname)

fprintf('End: mu=%f\n',...
    mu);

end

