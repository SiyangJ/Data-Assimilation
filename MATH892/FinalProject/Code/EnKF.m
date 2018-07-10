%% Summary of document

% Author: Siyang Jing
% Time: Mon Apr 16
% MATH892 Final Project
% Reference code: Amit Apte's ensemble kalman filter
% Reference paper: primarily Baek 2006

% Files to use together:
% L96: lorenz-96 model
% L96ME: lorenz-96 model with model error zeta and xi
% aux: for evaluation purposes

% Summary:
% This is version 3;
% For now, no error discovered;
% Version 1 and 2 can be regarded as incorrect;
% L96B1,2,3 are not used;
% They are not consistent with the paper I'm referring to;
% However, potentially interesting to try.

% TODO List:
% Need to reduce code redundancy!!!!
% Need to optimize the details of the algorithm.
% Need to inspect Type C model error, i.e. A,B both nonzero
    % inflmu 1e-3, 1e-4, are absolutely horrible.
    % inflmu 1 is way better.
    % Seems that my results for various inflmu is different from the paper
    % Need to solve other issues and set up the big test
% Need to function-ize this document to allow change of critical parameters.
% Which parameters to change?
% inflmu, according to the original paper;
% nens, which is too large here;
% nobs and deltaobs, to find stability of algorithm;
% A and B, in case of larger error?

% Interesting things to do:
% Fm is essentially not changed.
% We could change the forecast dynamical system completely, like
% x' = L(x) + zeta
% Something like this
% Intuitively, this could allow larger deltaobs?

%% Parameters allowed for change

F = 8;
A = 0.2*F;
B = 0.2*F;

nobs = 200;
deltaobs = 0.05;

ndim = 40;

dobs = 40;

nens = 40;

inflmu = 1;

%% Other parameters

rng(1);


sinus = sin(2*pi/ndim*(0:(ndim-1))');
zeta = A*sinus;
xi = B*sinus;

% True model
TM = @(t,x) L96ME(t,x,F,xi,zeta);
% Forecast models
TF = @(t,x) L96(t,x,F);

% Wrong!!!
%TFb1 = @(t,x) L96B1(t,x,F,ndim);
%TFb2 = @(t,x) L96B2(t,x,F,ndim);
%TFb3 = @(t,x) L96B3(t,x,F,ndim);

% Observation operators
H = [eye(dobs),zeros(dobs,ndim-dobs)];

Hb1 = [eye(dobs),zeros(dobs,2*ndim-dobs)];
Hb2 = [eye(dobs),zeros(dobs,ndim-dobs),eye(dobs),zeros(dobs,ndim-dobs)];
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

%% Plot a long trajectory

tlong = linspace(0,100,1e4);
% //TOCHANGE 
[~,xlong] = ode45(TM,tlong,xattr);
xlong = xlong';

%sigmaobs = (max(xlong(1,:))-min(xlong(1,:)))/50;
%sigmainit = sigmaobs*10;

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

xfm = zeros(ndim,nobs+1);
xam = zeros(ndim,nobs+1);
xfcov = zeros(ndim,ndim,nobs+1);
xacov = zeros(ndim,ndim,nobs+1);

xmf = zeros(ndim,nobs+1);
xcovf = zeros(ndim,ndim,nobs+1);

xfmb1 = zeros(2*ndim,nobs+1);
xamb1 = zeros(2*ndim,nobs+1);
xfcovb1 = zeros(2*ndim,2*ndim,nobs+1);
xacovb1 = zeros(2*ndim,2*ndim,nobs+1);

xfmb2 = zeros(2*ndim,nobs+1);
xamb2 = zeros(2*ndim,nobs+1);
xfcovb2 = zeros(2*ndim,2*ndim,nobs+1);
xacovb2 = zeros(2*ndim,2*ndim,nobs+1);

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

xamb2(1:ndim,1) = xam(:,1);
xacovb2(1:ndim,1:ndim,1) = xacov(:,:,1);
xacovb2(ndim+1:end,ndim+1:end,1) = sigmaxi^2*eye(ndim);

xamb3(1:ndim,1) = xam(:,1);
xacovb3(1:ndim,1:ndim,1) = xacov(:,:,1);
xacovb3(ndim+1:2*ndim,ndim+1:2*ndim,1) = sigmazeta^2*eye(ndim);
xacovb3(2*ndim+1:end,2*ndim+1:end,1) = sigmaxi^2*eye(ndim);

% Freerun ensemble
xensf = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (perfect model assumption)
xens = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Initial ensemble (bias model 1)
xensb1 = mvnrnd(xamb1(:,1),xacovb1(:,:,1),nens)';

% Initial ensemble (bias model 2)
xensb2 = mvnrnd(xamb2(:,1),xacovb2(:,:,1),nens)';

% Initial ensemble (bias model 3)
xensb3 = mvnrnd(xamb3(:,1),xacovb3(:,:,1),nens)';

%% Perform EnKF

for i=2:nobs+1
    
    % Forecast
    % Evolve the ensemble
    for j = 1:nens   
        
        [~,xf] = ode45(TF,tobs(i-1:i),xensf(:,j));
        xf = xf(end,:)';
        xensf(:,j) = xf;
        
        [~,xf] = ode45(TF,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xens(:,j) = xf; 
        
        % Bias Model 1
        [~,xf] = ode45(TF,tobs(i-1:i),xensb1(1:ndim,j));
        xf = xf(end,:)';
        xensb1(1:ndim,j) = xf + xensb1(ndim+1:end,j); 
        
        % Bias Model 2
        [~,xf] = ode45(TF,tobs(i-1:i),xensb2(1:ndim,j));
        xf = xf(end,:)';
        xensb2(1:ndim,j) = xf; 
        
        % Bias Model 3
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
    
    xfmb2(:,i) = mean(xensb2,2);
    xfcovb2(:,:,i) = cov(xensb2');
    
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
    pf = xfcovb2(:,:,i)+inflmu*trace(xfcovb2(:,:,i))/2/ndim*eye(2*ndim);
    pfht = pf*Hb2';   
    kgainb2 = pfht*(Hb2*pfht+Robsmat)^-1;    
    
    % Bias Model 3
    pf = xfcovb3(:,:,i)+inflmu*trace(xfcovb3(:,:,i))/3/ndim*eye(3*ndim);
    pfht = pf*Hb3';   
    kgainb3 = pfht*(Hb3*pfht+Robsmat)^-1;   
    
    % Update Ensemble members
    
    for j=1:nens     
        % For all models, use the same perturbed observations
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xens(:,j) = xens(:,j)+kgain*(pertob-H*xens(:,j));  
        xensb1(:,j) = xensb1(:,j)+kgainb1*(pertob-Hb1*xensb1(:,j));       
        xensb2(:,j) = xensb2(:,j)+kgainb2*(pertob-Hb2*xensb2(:,j));       
        xensb3(:,j) = xensb3(:,j)+kgainb3*(pertob-Hb3*xensb3(:,j));       
    end  
    
    % Calculate the stats of ensemble for xa    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
    
    xamb1(:,i) = mean(xensb1,2);
    xacovb1(:,:,i) = cov(xensb1');
    
    xamb2(:,i) = mean(xensb2,2);
    xacovb2(:,:,i) = cov(xensb2');
    
    xamb3(:,i) = mean(xensb3,2);
    xacovb3(:,:,i) = cov(xensb3');
    
end


%% Other stuff

for m = 0:3
figure
for i=20:23
    subplot(1,4,i-19)
    hold on
    plot(ttraj,truetraj(i,:),'b-','LineWidth',4)
    if m==0
        ami = xam(i,:);
        fmi = xfm(i,:);
        ame = reshape(xacov(i,i,:),[1,nobs+1]);
        fme = reshape(xfcov(i,i,:),[1,nobs+1]);
    elseif m==1
        ami = xamb1(i,:);
        fmi = xfmb1(i,:);
        ame = reshape(xacovb1(i,i,:),[1,nobs+1]);
        fme = reshape(xfcovb1(i,i,:),[1,nobs+1]);
    elseif m==2
        ami = xamb2(i,:)+xamb2(ndim+i,:);
        fmi = xfmb2(i,:)+xfmb2(ndim+i,:);
        ame = reshape(xacovb2(i,i,:),[1,nobs+1]);
        fme = reshape(xfcovb2(i,i,:),[1,nobs+1]);
    elseif m==3
        ami = xamb3(i,:)+xamb3(2*ndim+i,:);
        fmi = xfmb3(i,:)+xfmb3(2*ndim+i,:);
        ame = reshape(xacovb3(i,i,:),[1,nobs+1]);
        fme = reshape(xfcovb3(i,i,:),[1,nobs+1]);
    end
    errorbar(tobs,ami,ame,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    errorbar(tobs,fmi,fme,'r-s','LineWidth',2,'MarkerSize',8)
    plot(tobs,xmf(i,:),'g--*')
    hold off
end
end

