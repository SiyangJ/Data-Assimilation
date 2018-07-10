% Linear Shallow Water Model Error Client

close all;

T_i = 0;    T_f = 10;                                  % initial and final times
h = 0.01;                                              % computational time step
comp_time = T_i:h:T_f;                        % computational time domain
comp_len = length(comp_time);
ic = [0.4; 0.4; 1; 1];                                  % initial condition

% Generate the truth using RK4 solver
truth = zeros(4,comp_len);
truth(:,1) = ic;
for j=1:comp_len-1
    truth(:,j+1) = rk4('ens_lsw_ti_inertia',comp_time(j),comp_time(j+1),h,truth(:,j));
end

% Set parameters for observations
obs_int = 1;                                            % time interval between observations
obs_times = T_i:obs_int:T_f;                            % observation times
obs_len = length(obs_times);                            % number of observations
time_rat = obs_int/h;                                   % ratio of obs time interval to comp time interval
obs_var = 0.0001;                                         % observation error variance
R = obs_var*eye(2);
obs = zeros(2,obs_len);                               % matrix to hold observations
H = [1 0 0 0;
    0 1 0 0];                                           % observation operator

% Generate observations
for j=1:obs_len
    obs(:,j) = H*truth(:,time_rat*(j-1)+1)+normrnd(0,sqrt(obs_var),2,1);
end

%% Ensemble Kalaman Filter
ens_size = 100;                                        % number of ensemble members
ens_var = 0.01;                                      % initial ensemble variance
a_mean = zeros(4,obs_len);                        % tracks ensemble mean over time
f_mean = zeros(4,obs_len);
ensemble = repmat(truth(:,1),1,ens_size)+normrnd(0,sqrt(ens_var),4,ens_size);     % initial ensemble
for i=1:ens_size
    a_mean(:,1) = a_mean(:,1)+ensemble(:,i)/ens_size;           % compute ensemble mean
end

f_mean(:,1) = a_mean(:,1);

% Carry out EnKF
P_mats = zeros(4,4*obs_len);
K_mats = zeros(4,2*obs_len);
obs_error = zeros(1,obs_len);
xy_error = zeros(1,obs_len);
xy_error(1) = norm(a_mean(1:2,1)-truth(1:2,1));
inflation = 1;
fc_ensemble = zeros(4,ens_size);
an_ensemble = zeros(4,ens_size);

for j=1:obs_len-1
   
    % Forecast ensemble members forward in time
    for k=1:ens_size
        fc_ensemble(:,k) = rk4('ens_lsw_ti',obs_times(j),obs_times(j+1),h,ensemble(:,k));
    end
    %fore = rk4_full('ens_lsw_ti',obs_times(j),obs_times(j+1),h,a_mean(:,j));
    %plot(fore(1,:),fore(2,:),'g')
    for i=1:ens_size
        f_mean(:,j+1) = f_mean(:,j+1)+fc_ensemble(:,i)/ens_size;
    end
    
    % Compute sample forecast covariance matrix
    P = zeros(4,4);
    for i=1:ens_size
       P = P+(ensemble(:,i)-f_mean(:,j+1))*transpose(ensemble(:,i)-f_mean(:,j+1))/(ens_size-1);
    end
    
    P = inflation*P;
    
    % Compute Kalman gain
    K = P*H'*pinv(H*P*H'+R);
    
    % Update ensemble members
    for i=1:ens_size
        an_ensemble(:,i) = fc_ensemble(:,i)+K*(obs(:,j+1)+normrnd(0,sqrt(obs_var),2,1)-H*fc_ensemble(:,i));
    end
    
    % Compute ensemble mean
    for i=1:ens_size
        a_mean(:,j+1) = a_mean(:,j+1)+an_ensemble(:,i)/ens_size;
    end
    
    ensemble = an_ensemble;
    
    obs_error(j+1) = sqrt((obs(1,j+1)-truth(1,time_rat*j+1))^2+(obs(2,j+1)-truth(2,time_rat*j+1))^2)/sqrt(obs_var);
    xy_error(j+1) = sqrt(((a_mean(1,j+1)-truth(1,time_rat*j+1))^2+(a_mean(2,j+1)-truth(2,time_rat*j+1))^2));
    
end


%% Plot results

figure
hold on
plot(truth(1,:),truth(2,:))
scatter(obs(1,:),obs(2,:),'r')
scatter(f_mean(1,2:end),f_mean(2,2:end),'g')
scatter(a_mean(1,1:end),a_mean(1,1:end),'k')
xlabel('x')
ylabel('y')
title('Perfect Model')
legend('Truth','Observations','Forecast','Analysis')

figure
hold on
plot(obs_times,xy_error)
avg_error = cumsum(xy_error)/length(obs_times);
plot(obs_times,avg_error,'r')
legend('EnKF Error','Avg Observational Error')
xlabel('Time')
ylabel('RMS Error')
%title('Order 1 Model Error, Frequent Observations')

