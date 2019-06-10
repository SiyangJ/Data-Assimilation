%%

clc; close all; clear;

%% Parameters allowed for change

FILTER_TYPE = 3;

F = 8;

nobs = 200;
deltaobs = 0.05;

ndim = 40;

dobs = 20;

nens = 20;

inflmu = 0.05;

%% Other parameters

rng(1000);

% True model
TM = @(t,x) L96_10_30(t,x,F);

% Forecast models
TF = @(t,x) L96_10_30(t,x,F);

% Note: They could be different.
% That was something I was working on.
% But for now we just take the same.

% Observation operators
% H = [eye(dobs),zeros(dobs,ndim-dobs)];
% h = @(x) HML40_L_10_30(x);
h = @(x) obs_40_20(x);

h_jacobian = @(x) obs_40_20_jacobian(x);

% stupid inverse case
% h = @(x) 

% Seems like nonlinear not working, try linear
% h = @(x) x(1:20,:);
% H = [eye(dobs),zeros(dobs,ndim-dobs)];

% Magic number that works very well
sigmaobs = 0.09;

% Recommended choice by Baek 2006
sigmainit = 1.3;

%% Find an initial condition on the attractor

xrand = rand(ndim,1);
ttrans = linspace(0,100,1001);
[~,xtrans] = ode45(TM,ttrans,xrand);

xtrans = xtrans';
xattr = xtrans(:,end);

%% Plot a long trajectory
%  This part is not actually used by the rest part of the code.
%  Amit suggested taking initial spread and observation covariance from the long trajectory, as coded below.

tlong = linspace(0,100,1e4);
[~,xlong] = ode45(TM,tlong,xattr);
xlong = xlong';

% This method was used by Amit
% sigmaobs = (max(xlong(1,:))-min(xlong(1,:)))/50;
% sigmainit = sigmaobs*10;

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
% Change from linear to full
% trueobs = H*truestate;
trueobs = h(truestate);
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

% Initial mean and covariance
% Need to test initial mean
xam(:,1) = xattr + sigmainit*randn(ndim,1);
xacov(:,:,1) = sigmainit^2*eye(ndim);

xmf(:,1) = xam(:,1);
xcovf(:,:,1) = xacov(:,:,1);

% Freerun initial ensemble
xensf = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

% Model initial ensemble
xens = mvnrnd(xam(:,1),xacov(:,:,1),nens)';

%% Perform EnKF

%for i=2:nobs+1
    
    % Forecast
    % Evolve the ensemble
    for j = 1:nens   
        
        [~,xf] = ode45(TF,tobs(i-1:i),xensf(:,j));
        xf = xf(end,:)';
        xensf(:,j) = xf;
        
        [~,xf] = ode45(TF,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xf = xf + randn(40,1);
        xens(:,j) = xf; 
        
    end
    
    % Calculate the stats of ensemble for xf
   
    xmf(:,i) = mean(xensf,2);
    xcovf(:,:,i) = cov(xens');
    
    xfm(:,i) = mean(xens,2);
    xfcov(:,:,i) = cov(xens');
    
    yens = h(xens);
    yfm  = mean(yens,2);
    
    % Update
    
    % Kalman gain
    % Fully estimated by sample covariances
    
    if FILTER_TYPE == 1
    
    
        % Old methods, not fully estimated.
        % Enhanced variance inflation
        % Ott et al. 2004
        pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
        pfht = pf*H';   
        kgain = pfht*pinv(H*pfht+Robsmat);  

    elseif FILTER_TYPE == 2
    
        Ex  = xens-xfm(:,i);
        Ey  = yens-yfm;
        Bxy = Ex*Ey'/(nens-1);
        Byy = cov(yens') + Robsmat;
        kgain = Bxy*pinv(Byy);
    
    elseif FILTER_TYPE == 3
        
        H = h_jacobian(xfm(:,i));
        pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
        pfht = pf*H';   
        kgain = pfht*pinv(H*pfht+Robsmat);  
        
    end
    
    % Update Ensemble members
    
    for j=1:nens     
        % For all models, use the same perturbed observations
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xens(:,j) = xens(:,j)+kgain*(pertob-h(xens(:,j)));       
    end  
    
    % Calculate the stats of ensemble for xa    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
    
    figure

    subplot(2,1,1)
    plot(1:10,truestate(1:10,i),'b-','LineWidth',6)
    hold on
    plot(1:10,xam(1:10,i),'r-','LineWidth',3)
    plot(1:10,xfm(1:10,i),'g-','LineWidth',1)
    
    hold off
    
    subplot(2,1,2)
    plot(11:40,truestate(11:40,i),'b-','LineWidth',6)
    hold on
    plot(11:40,xam(11:40,i),'r-','LineWidth',3)
    plot(11:40,xfm(11:40,i),'g-','LineWidth',1)
    
    hold off
    
%end
