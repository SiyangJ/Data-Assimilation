%%

clc; close all; clear;

%% Parameters allowed for change

INSPECT_ENSEMBLE_EVOLUTION = false;

FILTER_TYPE = 2;

nobs = 600;
deltaobs = 0.05;

ndim = 2;

dobs = 5;

nens = 20;

inflmu = 0.5;

%% Other parameters

rng(1000);

% True model
TM = @(t,x) EW09M(t,x,60,deltaobs);

% Forecast models
TF = @(t,x) EW09M(t,x,60,deltaobs);

% Note: They could be different.
% That was something I was working on.
% But for now we just take the same.

% Observation operators
% HM = [eye(dobs),zeros(dobs,ndim-dobs)];
% h = @(x) HML40_L_10_30(x);

% HM = zeros(20,40);
% for i = 1:20
%     HM(i,2*i-1:2*i) = [1,1];
% end
% 
% h = @(x) HM * x;

% HM = zeros(20,40);
% for i = 1:20
%     
%     HM(i,2*i-1) = 1;
% end

% HM = [eye(20),zeros(20,20)];
% 
% h = @(x) HM * x;
% 
% Hf = HM;
% hf = h;

% h = @(x) new_obs_40_20(x);
% hf = @(x) new_obs_40_20(x);

% h = @(x) EW09_Csat(x);
% 
% hf = @(x) EW09_Ci(x);

h = @(x) EW09_obs(x);

hf = @(x) EW09_obs(x);% + [1;0.1;1;100;0.01].*randn(dobs,size(x,2));


% h = @(x) obs_40_20(x);

% h_jacobian = @(x) obs_40_20_jacobian(x);

% h_jacobian = @(x) new_obs_40_20_jacobian(x);

% stupid inverse case
% h = @(x) 

% Seems like nonlinear not working, try linear
% h = @(x) x(1:20,:);
% HM = [eye(dobs),zeros(dobs,ndim-dobs)];

sigmaobs = [1;0.1;1;100;0.01];

% sigmaobs = [0.05];

% Recommended choice by Baek 2006
sigmainit = [3.2;0.1];

%% Find an initial condition on the attractor

xattr = [-200;0.7];

%% Plot a long trajectory
%  This part is not actually used by the rest part of the code.
%  Amit suggested taking initial spread and observation covariance from the long trajectory, as coded below.

%% Generate truth 

tend = nobs * deltaobs;
tobs = linspace(0,tend,nobs+1);

truestate = zeros(ndim,nobs+1);
truestate(:,1) = xattr;

for i=2:nobs+1
    
    truestate(:,i) = TM(deltaobs*(i-1),truestate(:,i-1));
    
end


%% Generate obervations

% Assume diagonal covariance
% Change from linear to full
% trueobs = HM*truestate;
trueobs = h(truestate);
robs = sigmaobs.^2;
Robsmat = diag(robs);
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
xam(:,1) = xattr + sigmainit.^2.*randn(ndim,1);
xacov(:,:,1) = sigmainit.^2.*eye(ndim);

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
        xlim([-200,100])
        ylim([0,0.9])
        pause();
    end
    
    % Forecast
    % Evolve the ensemble
    for j = 1:nens   
        
        xf = TF(tobs(i-1),xensf(:,j));
        xensf(:,j) = xf;
        
        xf = TF(tobs(i-1),xens(:,j));
        xens(:,j) = xf; 
        
    end
    
            
    if INSPECT_ENSEMBLE_EVOLUTION
        scatter(xens(1,:),xens(2,:));
        pause();
    end
    
    for j = 1:nens   
        
        xf = xensf(:,j);
        xf = xf + sigmainit.^2.*randn(ndim,1);
        xensf(:,j) = xf;
        
        xf = xens(:,j);
        xf = xf + sigmainit.^2.*randn(ndim,1);
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
        
        HM = h_jacobian(xfm(:,i));
        pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
        pfht = pf*HM';   
        kgain = pfht*pinv(HM*pfht+Robsmat);  
        
    end
    
    % Update Ensemble members
    
    for j=1:nens     
        % For all models, use the same perturbed observations
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
