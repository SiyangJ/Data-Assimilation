clear
close all

%% this section is for 'critical' parameters, that tend to be changed often.

Npf = 1000;    % number of particles

initial_error = 0.00; %std of initial error; scalar, or vector of length mdim
particle_spread = 0.1; %std of initial particle spread; set >= initial_error for 'good' prior

obs_step = .3; %time between observations

obs_var = 0.1; % observation variance

dt = 0.001; %model time step

resamp_thresh = 0.5; %threshold for sample diversity to trigger resampling
wiggle = 1e-1; %standard deviation of noise to add on resampling


%% Model/truth parameters

T = 14*obs_step; %final time

tvals = 0 : obs_step : T; %analysis times
tdim = length(tvals);

%%%%% stochastic bimodal model
% model = 'bimodal';
% truth=model;
% sigma=0.4; %standard deviation of model noise
% mdim=2; %can be any integer; model will then consist of mdim independent simulations
% tIC=zeros(mdim,1); %x=0 is an unstable equilibrium for the deterministic part
%
% H=eye(mdim); %observe every variable

% %%%%%lorenz63 model
model = 'lorenz63'; %set the file containing the model dyn. sys.
truth = 'lorenz63'; %could have different truth/models
mdim=3; %model dimension (fixed)
sigma=0.1; %standard deviation of model noise
tIC = rand(3,1); %true initial condition
% H = eye(mdim);  % observe every variable
 H = [1 1 1]/3; %observe the mean
% H = [1 0 0;
%      0 1 0]; %observe first two variables only

% %%%%%%%%%%%simple square model
% model='simplexsq';
% truth=model;
% mdim=4;
% sigma=0.1; %standard deviation of model noise
% tIC=rand(mdim,1);%(10,1);
% H=eye(mdim);

%%%%%%%%lorenz96
% model='lorenz96';
% truth='lorenz96';
% mdim=40; %model dimension (variable)
% tIC = 1*rand(mdim,1);
% H=eye(mdim);  % observation operator
% H=zeros(mdim);
% H(1,1)=1;
% sigma=0.1; %standard deviation of model noise


mIC = tIC + initial_error.*randn(mdim,1); %model initial guess

%% data parameters

[obsdim,~] = size(H); %dimension of observations
R = obs_var*eye(obsdim); %assuming observation errors independent


%% Generate ensemble

particle = repmat(mIC,1,Npf)+...
    normrnd(zeros(mdim,Npf),particle_spread);       % saves the 'current' particle ensemble, for calculations
size(particle)
phist = zeros(mdim,Npf,tdim); %history of particle values
phist(:,:,1) = particle;

W=ones(Npf,1)/Npf; %initial particle weights
Whist = zeros(Npf,tdim); %history of particle weights
Whist(:,1)= W;


%% Generate truth and observations
%deterministic model
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,xt] = ode45(truth,tvals,tIC,options); %'DA phrasing': xt for truth
xt = xt'; %order by x then t, to match ensemble ordering by x then member

%stochastic model
% tdt = 1e-4; %time step to use when generating the 'truth'
% xt = EulerM_hist(truth,tvals,tdt,tIC,sigma);

obshist = H*(xt + sqrt(obs_var)*randn(mdim,tdim));


%% Carry out PF

resampcount=0; %count how many times the PF resampled

for tau=1:tdim-1
    
    % update particles
    
    %deterministic model
%     [particle] = rk4(model,[tvals(tau) tvals(tau+1)],dt,particle);
    
    %stochastic model
     particle = EulerM(model,tvals(tau:tau+1),dt,particle,sigma);
    
    phist(:,:,tau+1) = particle;
    
    % get observation
    obs = obshist(:,tau+1); %we never use the observation at the initial time; by convention the first observation is at t_1, not the initial time.
    
    %update weights - Gaussian assumption for error statistics
    innov = obs-H*particle; %innovation between observation and 'the predicted observation' H*x_f
    %note - if weights are updated using W=W*exp(-innov), as written in
    %textbooks, then W=0 everywhere is a common problem. The calculations
    %below avoid that problem.
    Wtmp =  -0.5*innov'/R*innov;    Wmax=max(Wtmp);
    Wtmp  = Wtmp-Wmax;
    W = W.*exp(Wtmp'); W=W/sum(W);
    Whist(:,tau+1) = W;
    
    %test for resampling
    if   1/sum(W.^2)/Npf < resamp_thresh %test for resampling
        resampcount=resampcount+1;
        sampIndex = resampleMultinomial(W);
        particle = particle(:,sampIndex)+wiggle*randn(mdim,Npf);
        W = ones(Npf,1)/Npf;
        Whist(:,tau+1)=W;
    end
    
end


%% Plot results
[a,b]=findIntegerFactors(tdim);

% figure
% for ind=1:tdim
%     subplot(a,b,ind)
%     hold on
%     title(['t= ' num2str(tvals(ind))])
%     %axis([-2 2 -2 2])
%     scatter(phist(1,:,ind),phist(2,:,ind),.5,'k') %plot all particles at uniform size
%     plot(xt(1,ind),xt(2,ind),'kx','LineWidth',2,'MarkerSize',10) %true state value
%     %plot(obshist(1,ind),obshist(2,ind),'rx','LineWidth',2,'MarkerSize',10) %observed state value
% end
% subtitle('Unweighted particle distribution')

figure
for ind=1:tdim
    subplot(a,b,ind)
    hold on
    title(['t= ' num2str(tvals(ind))])
    %axis([-2 2 -2 2])
    scatter(phist(1,:,ind),phist(2,:,ind),.01+5*Whist(:,ind)/max(Whist(:,ind)),'b') %plot particles with size given by weight
    plot(xt(1,ind),xt(2,ind),'kx','LineWidth',2,'MarkerSize',10) %true state value
    %plot(obshist(1,ind),obshist(2,ind),'rx','LineWidth',2,'MarkerSize',10) %observed state value
end
subtitle(sprintf('%s,Weighted particle distribution',model))







