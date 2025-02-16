%%

clc; close all; clear;

%% Parameters allowed for change

INSPECT_ENSEMBLE_EVOLUTION = true;

FILTER_TYPE = 1;

F = 8;

nobs = 200;
deltaobs = 0.05;

ndim = 40;

dobs = 10;

nens = 20;

inflmu = 0.5;

%% Other parameters

rng(1000);

% True model
TM = @(t,x) L96(t,x,F);

% Forecast models
TF = @(t,x) L96(t,x,F);

% Note: They could be different.
% That was something I was working on.
% But for now we just take the same.

% Observation operators
% H = [eye(dobs),zeros(dobs,ndim-dobs)];
% h = @(x) HML40_L_10_30(x);

% H = zeros(20,40);
% for i = 1:20
%     H(i,2*i-1:2*i) = [1,1];
% end
% 
% h = @(x) H * x;

% H = zeros(20,40);
% for i = 1:20
%     
%     H(i,2*i-1) = 1;
% end

% H = [eye(20),zeros(20,20)];
% 
% h = @(x) H * x;
% 
% Hf = H;
% hf = h;

% h = @(x) new_obs_40_20(x);
% hf = @(x) new_obs_40_20(x);

% h = @(x) new_stupid_inverse(x);
% 
% Hf = [eye(10),zeros(10,30)];
% hf = @(x) Hf * x;
% 
% H = zeros(8,40);
% for i = 1:8
%     H(i,(i-1)*5+1:i*5) = 1/5;
% end
% 
% Hf = H;

% H = [eye(8),zeros(8,32)];

% H = zeros(10,40);
% for i = 1:10
%     H(i,(i-1)*4+1:i*4) = [5,1,10,2]/18;
% end

% H = zeros(2,40);
% H(1,1:5) = [5,1,10,2,20];
% H(2,end-5:end) = [2,4,9,20,1];

% H = [[5,1,10,2,20],zeros(1,35)];

% H = zeros(13,39);
% for i = 1:13
%     H(i,(i-1)*3+1:i*3) = [5,1,10]/16;
% end

% H = zeros(2,ndim);
% H(1,1:3) = [5,1,10]/16;
% H(2,4:6) = [2,4,20]/26;

% Hf = H;

% H = eye(40);
% 
% H = zeros(5,40);
% 
% for i = 1:5
%     H(i,2*i-1:2*i) = 1;
% end

H = [eye(10),zeros(10,30)];
Hf = H;

h = @(x) new_stupid_inverse(x);
hf = @(x) Hf * x;

% h = @(x) obs_40_20(x);

% h_jacobian = @(x) obs_40_20_jacobian(x);

% h_jacobian = @(x) new_obs_40_20_jacobian(x);

% stupid inverse case
% h = @(x) 

% Seems like nonlinear not working, try linear
% h = @(x) x(1:20,:);
% H = [eye(dobs),zeros(dobs,ndim-dobs)];

sigmaobs = 0.09;

% Recommended choice by Baek 2006
sigmainit = 1;

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
if INSPECT_ENSEMBLE_EVOLUTION
    figure
end
for i=2:nobs+1
    
    fprintf('Started i=%d',i);
    
    if INSPECT_ENSEMBLE_EVOLUTION
        clf
        scatter(xens(1,:),xens(2,:));
        hold on
        scatter(truestate(1,i),truestate(2,i),'Marker','*')
        xlim([-15,15])
        ylim([-15,15])
        pause();
    end
    
    % Forecast
    % Evolve the ensemble
    for j = 1:nens   
        
        [~,xf] = ode45(TF,tobs(i-1:i),xensf(:,j));
        xf = xf(end,:)';
        xensf(:,j) = xf;
        
        [~,xf] = ode45(TF,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xf = xf + sigmainit*randn(ndim,1);
        xens(:,j) = xf; 
        
    end
    
    % Calculate the stats of ensemble for xf
   
    xmf(:,i) = mean(xensf,2);
    xcovf(:,:,i) = cov(xens');
    
    xfm(:,i) = mean(xens,2);
    xfcov(:,:,i) = cov(xens');
    
        
    if INSPECT_ENSEMBLE_EVOLUTION
        scatter(xens(1,:),xens(2,:));
        scatter(xfm(1,i),xfm(2,i));
        pause();
    end
    
    yens = hf(xens);
    yfm  = mean(yens,2);
    
    % Update
    
    % Kalman gain
    % Fully estimated by sample covariances
    if FILTER_TYPE < 4
    if FILTER_TYPE == 1
    
    
        % Old methods, not fully estimated.
        % Enhanced variance inflation
        % Ott et al. 2004
        pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
        pfht = pf*Hf';   
        kgain = pfht*pinv(Hf*pfht+Robsmat);  

    elseif FILTER_TYPE == 2
    
%         IC = inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
%         
%         for j = 1:nens
%             
%             xens(:,j) = mvnrnd(xens(:,j), IC);
%         end

%         pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
%         
%         xens_new = mvnrnd(xfm(:,i), pf, nens)';

        xens_new = xens;
        if INSPECT_ENSEMBLE_EVOLUTION
        scatter(xens_new(1,:),xens_new(2,:));
        
        pause();
        end
        yens = hf(xens_new);
        yfm = mean(yens,2);
        
        Ex  = xens_new-xfm(:,i);
        % Ex = Ex * (1+inflmu);
        
        Ey  = yens-yfm;
        Bxy = Ex*Ey'/(nens-1);
        Byy = cov(yens');
        Byy = Byy + inflmu*trace(Byy)/dobs*eye(dobs) + Robsmat;
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
        % disp(yobs(:,i) - hf(xens(:,j)));
        pertob = mvnrnd(yobs(:,i),Robsmat,1)';
        xens(:,j) = xens(:,j)+kgain*(pertob-hf(xens(:,j)));       
    end
    
    elseif FILTER_TYPE == 4
        
        % ETKF
        
        Ex  = xens-xfm(:,i);
        Ey  = yens-yfm;
        
        C = Ey'*Robsmat^-1;   
        
        % inflation not useful ... ?
        Ptil = pinv((nens-1) * eye(nens) / (1+inflmu) + C*Ey);        
        Wa = sqrtm((nens-1)*Ptil);
        
        wa = Ptil * C * (yobs(:,i) - yfm);
        
        % Wa is for perturbation/anomaly in the ensemble
        % wa is for correction from forecast to analysis
        
        W = Wa + wa;
        
        xens = Ex*W + xfm(:,i);
        
    elseif FILTER_TYPE == 5
        
        % ETKF 2
        Linv = chol(Robsmat)^-1;
        S1 = (xens-xfm(:,i)) / sqrt(nens-1);
        F = (yens-yfm)' * Linv' / sqrt(nens-1);
        T = F' * F + eye(dobs);
        
        K = S1 * F * T^-1;
        
        xam = xfm(:,i) + K*(Linv * (yobs(:,i) - hf(xfm(:,i))));
        
        [V,L] = eig(F*F');
        
        S2 = S1 * V * sqrtm(L+eye(nens))^-1 * V';
        
        xens = xam + S2 / sqrt(nens-1);
    
    end
    % Calculate the stats of ensemble for xa    
    xam(:,i) = mean(xens,2);
    xacov(:,:,i) = cov(xens');
    
    if INSPECT_ENSEMBLE_EVOLUTION
        scatter(xens(1,:),xens(2,:));
        scatter(xam(1,i),xam(2,i))
        hold off
        pause()
    end
    
    fprintf('Finished i=%d\n',i);
    
end
