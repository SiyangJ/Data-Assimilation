function BM = EnKFfun_LSW_v1(mu,epsilon,fname)

DEBUG__ = true;

%% Parameters allowed for change

fprintf('Start: mu=%f\n',...
    mu);

nobs = 500;
deltaobs = 0.05;

ndim = 4;

dobs = 2;

% Seems a reasonable choice
nens = 20;

inflmu = mu;

%% Other parameters

rng(1);

% True model
TM = @(t,x) lsw_ti_inertia(t,x);
% Forecast models
% Perfect
TF1 = @(t,x) lsw_ti_inertia(t,x);
% Linear
TF2 = @(t,x) lsw_ti(t,x);
% Order1
TF3 = @(t,x) lsw_ti_order1(t,x);

% Observation operators
% For Lagrangian DA, just the driftors
H = [1 0 0 0;...
     0 1 0 0];

sigmaobs = 0.09;

sigmainit = 1.3;

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

% Perfect
xfm = zeros(ndim,nobs+1);
xam = zeros(ndim,nobs+1);
xfcov = zeros(ndim,ndim,nobs+1);
xacov = zeros(ndim,ndim,nobs+1);

% Free run
xmf = zeros(ndim,nobs+1);
xcovf = zeros(ndim,ndim,nobs+1);

% Linear
xfmb1 = zeros(ndim,nobs+1);
xamb1 = zeros(ndim,nobs+1);
xfcovb1 = zeros(ndim,ndim,nobs+1);
xacovb1 = zeros(ndim,ndim,nobs+1);

% Order1
xfmb2 = zeros(ndim,nobs+1);
xamb2 = zeros(ndim,nobs+1);
xfcovb2 = zeros(ndim,ndim,nobs+1);
xacovb2 = zeros(ndim,ndim,nobs+1);

% Initial mean and covariance
% Need to test initial mean
xam(:,1) = xattr + sigmainit*randn(ndim,1);
xacov(:,:,1) = sigmainit^2*eye(ndim);

xmf(:,1) = xam(:,1);
xcovf(:,:,1) = xacov(:,:,1);

% For bias models
xamb1(:,1) = xam(:,1);
xacovb1(:,:,1) = xacov(:,:,1);

xamb2(:,1) = xam(:,1);
xacovb2(:,:,1) = xacov(:,:,1);

% Freerun ensemble
xensf  = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (perfect model)
xens   = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble bias models
xensb1 = mvnrnd(xamb1(:,1),xacovb1(:,:,1),nens)';
xensb2 = mvnrnd(xamb2(:,1),xacovb2(:,:,1),nens)';

%% Perform EnKF

for i=2:nobs+1
    
    % Forecast
    % Evolve the ensemble
    
    % Free run and Perfect Model
    for j = 1:nens   
        % Use perfect for free run
        [~,xf] = ode45(TF1,tobs(i-1:i),xensf(:,j));
        xf = xf(end,:)';
        xensf(:,j) = xf;
        
        % Perfect model
        [~,xf] = ode45(TF1,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xens(:,j) = xf; 
        
        % Linear
        [~,xf] = ode45(TF2,tobs(i-1:i),xensb1(:,j));
        xf = xf(end,:)';
        xensb1(:,j) = xf; 
        
        % Order 1 Inertia
        [~,xf] = ode45(TF3,tobs(i-1:i),xensb2(:,j));
        xf = xf(end,:)';
        xensb2(:,j) = xf; 
        
    end
    
    % Calculate the stats of ensemble for xf
   
    xmf(:,i) = mean(xensf,2);
    xcovf(:,:,i) = cov(xens');
    
    xfm(:,i) = mean(xens,2);
    xfcov(:,:,i) = cov(xens');
    
    xfmb1(:,i) = mean(xensb1,2);
    xfcovb1(:,:,i) = cov(xensb1');
    
    xfmb2(:,i) = mean(xensb2,2);
    xfcovb2(:,:,i) = cov(xensb2');
    
    % Update
    
    % Kalman gain
    
    % Perfect Model   
    % Enhanced variance inflation
    % Ott et al. 2004
    pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
    pfht = pf*H';   
    kgain = pfht*(H*pfht+Robsmat)^-1;   
    
    % Linear  
    pf = xfcovb1(:,:,i)+inflmu*trace(xfcovb1(:,:,i))/ndim*eye(ndim);
    pfht = pf*H';   
    kgainb1 = pfht*(H*pfht+Robsmat)^-1;   
    
    % Order 1
    pf = xfcovb2(:,:,i)+inflmu*trace(xfcovb2(:,:,i))/ndim*eye(ndim);
    pfht = pf*H';   
    kgainb2 = pfht*(H*pfht+Robsmat)^-1;     
    
    % Update Ensemble members
    
    for j=1:nens     
        % Perturb the observations
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        % Apply Kalman Filter
        xens(:,j)   = xens(:,j)  +kgain  *(pertob-H*xens(:,j));   
        xensb1(:,j) = xensb1(:,j)+kgainb1*(pertob-H*xensb1(:,j));            
        xensb2(:,j) = xensb2(:,j)+kgainb2*(pertob-H*xensb2(:,j));
    end  
        
    % Calculate the stats of ensemble for xa    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
    
    xamb1(:,i) = mean(xensb1,2);
    xacovb1(:,:,i) = cov(xensb1');
    
    xamb2(:,i) = mean(xensb2,2);
    xacovb2(:,:,i) = cov(xensb2');
    
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
sumrmsab2 = 0;

tspan = floor(nobs/2)+1:nobs+1;
T = length(tspan);

for t = tspan
    
    sumrmsap = sumrmsap + norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsab1 = sumrmsab1 + norm(xamb1(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsab2 = sumrmsab2 + norm(xamb2(:,t)-truestate(:,t))/sqrt(ndim);
    
end

eap = sumrmsap/T;
eab1 = sumrmsab1/T;
eab2 = sumrmsab2/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

sumrmsp = 0;
sumrmsb1 = 0;
sumrmsb2 = 0;

for t = tspan
    
    sumrmsp = sumrmsp + norm(xfm(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb1 = sumrmsb1 + norm(xfmb1(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb2 = sumrmsb2 + norm(xfmb2(:,t)-truestate(:,t))/sqrt(ndim);
    
end

ep = sumrmsp/T;
eb1 = sumrmsb1/T;
eb2 = sumrmsb2/T;

BM = {inflmu,epsilon,deltaobs,nens,...
     eap,eab1,eab2,...
     ep,eb1,eb2};

save(fname)

fprintf('End: mu=%f\n',...
    mu);

end

