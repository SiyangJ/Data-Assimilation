% Model error client

close all;

T_i = 0;    T_f = 50;
h = 0.01;
comp_time = T_i:h:T_f;
comp_len = length(comp_time);

obs_int = 1;
obs_times = T_i:obs_int:T_f;
obs_len = length(obs_times);
time_rat = obs_int/h;
obs_var = 0.01;
R = obs_var*eye(2);
H = [1 0 0 0;
    0 1 0 0];

avg_perf_error = zeros(1,obs_len);
avg_noninert_error = zeros(1,obs_len);
avg_order1_error = zeros(1,obs_len);
avg_obs_error = zeros(1,obs_len);

num_trials = 25;
ens_size = 50;

for i=1:num_trials
    
    truth = zeros(4,comp_len);
    %ic = [unifrnd(0,0.5); unifrnd(0,0.5); 1; 1];
    ic = [0.4; 0.25; 1; 1];
    truth(:,1) = ic;
    for j=1:comp_len-1
        truth(:,j+1) = rk4('ens_lsw_ti_inertia',comp_time(j),comp_time(j+1),h,truth(:,j));
    end
    
    obs = zeros(2,obs_len);
    for j=1:obs_len
        obs(:,j) = H*truth(:,time_rat*(j-1)+1)+normrnd(0,sqrt(obs_var),2,1);
    end
    
    obs_error = zeros(1,obs_len);
    for j=1:obs_len
        obs_error(j) = sqrt((obs(1,j)-truth(1,time_rat*(j-1)+1))^2+(obs(2,j)-truth(2,time_rat*(j-1)+1))^2);
    end
    avg_obs_error = avg_obs_error+obs_error/num_trials;
    
    perf_ensemble = zeros(4*ens_size,obs_len);                  % ensembles
    noninert_ensemble = zeros(4*ens_size,obs_len);
    order1_ensemble = zeros(4*ens_size,obs_len);
    ens_var = 0.1;                                              % initial ensemble variance
    
    perf_fmean = zeros(4,obs_len);                              % forecast and analysis means
    perf_amean = zeros(4,obs_len);
    noninert_fmean = zeros(4,obs_len);
    noninert_amean = zeros(4,obs_len);
    order1_fmean = zeros(4,obs_len);
    order1_amean = zeros(4,obs_len);
    
    perturb = [normrnd(0,sqrt(obs_var)); normrnd(0,sqrt(obs_var)); 0; 0];
    perf_ensemble(:,1) = repmat(truth(:,1)+perturb,ens_size,1)+normrnd(0,sqrt(ens_var),4*ens_size,1);
    noninert_ensemble(:,1) = repmat(truth(:,1)+perturb,ens_size,1)+normrnd(0,sqrt(ens_var),4*ens_size,1);
    order1_ensemble(:,1) = repmat(truth(:,1)+perturb,ens_size,1)+normrnd(0,sqrt(ens_var),4*ens_size,1);
    
    for k=1:ens_size
        perf_amean(:,1) = perf_amean(:,1)+perf_ensemble(4*k-3:4*k,1)/ens_size;
        noninert_amean(:,1) = noninert_amean(:,1)+noninert_ensemble(4*k-3:4*k,1)/ens_size;
        order1_amean(:,1) = order1_amean(:,1)+order1_ensemble(4*k-3:4*k,1)/ens_size;
    end
    
    avg_perf_error(1) = avg_perf_error(1)+norm(perf_amean(1:2,1)-truth(1:2,1))/num_trials;
    avg_noninert_error(1) = avg_noninert_error(1)+norm(noninert_amean(1:2,1)-truth(1:2,1))/num_trials;
    avg_order1_error(1) = avg_order1_error(1)+norm(order1_amean(1:2,1)-truth(1:2,1))/num_trials;
    
    for j=1:obs_len-1
        
        perf_ensemble(:,j+1) = rk4('ens_lsw_ti_inertia',obs_times(j),obs_times(j+1),h,perf_ensemble(:,j));
        noninert_ensemble(:,j+1) = rk4('ens_lsw_ti',obs_times(j),obs_times(j+1),h,noninert_ensemble(:,j));
        order1_ensemble(:,j+1) = rk4('ens_lsw_ti_order1',obs_times(j),obs_times(j+1),h,order1_ensemble(:,j));
        
        for k=1:ens_size
            perf_fmean(:,j+1) = perf_fmean(:,j+1)+perf_ensemble(4*k-3:4*k,j+1)/ens_size;
            noninert_fmean(:,j+1) = noninert_fmean(:,j+1)+noninert_ensemble(4*k-3:4*k,j+1)/ens_size;
            order1_fmean(:,j+1) = order1_fmean(:,j+1)+order1_ensemble(4*k-3:4*k,j+1)/ens_size;
        end
        
        P_perf = zeros(4,4);    P_noninert = zeros(4,4);    P_order1 = zeros(4,4);
        
        % forecast covariances
        for k=1:ens_size
            P_perf = P_perf+(perf_ensemble(4*k-3:4*k,j+1)-perf_fmean(:,j+1))...
                *transpose(perf_ensemble(4*k-3:4*k,j+1)-perf_fmean(:,j+1))/(ens_size-1);
            P_noninert = P_noninert+(noninert_ensemble(4*k-3:4*k,j+1)-noninert_fmean(:,j+1))...
                *transpose(noninert_ensemble(4*k-3:4*k,j+1)-noninert_fmean(:,j+1))/(ens_size-1);
            P_order1 = P_order1+(order1_ensemble(4*k-3:4*k,j+1)-order1_fmean(:,j+1))...
                *transpose(order1_ensemble(4*k-3:4*k,j+1)-order1_fmean(:,j+1))/(ens_size-1);
        end
        
        K_perf = P_perf*H'/(H*P_perf*H'+R);
        K_noninert = P_noninert*H'/(H*P_noninert*H'+R);
        K_order1 = P_order1*H'/(H*P_order1*H'+R);
        
        for k=1:ens_size
            perf_ensemble(4*k-3:4*k,j+1) = perf_ensemble(4*k-3:4*k,j+1)+...
                K_perf*(obs(:,j+1)+normrnd(0,sqrt(obs_var),2,1)-H*perf_ensemble(4*k-3:4*k,j+1));
            noninert_ensemble(4*k-3:4*k,j+1) = noninert_ensemble(4*k-3:4*k,j+1)+...
                K_noninert*(obs(:,j+1)+normrnd(0,sqrt(obs_var),2,1)-H*noninert_ensemble(4*k-3:4*k,j+1));
            order1_ensemble(4*k-3:4*k,j+1) = order1_ensemble(4*k-3:4*k,j+1)+...
                K_order1*(obs(:,j+1)+normrnd(0,sqrt(obs_var),2,1)-H*order1_ensemble(4*k-3:4*k,j+1));
        end
        
        for k=1:ens_size
            perf_amean(:,j+1) = perf_amean(:,j+1)+perf_ensemble(4*k-3:4*k,j+1)/ens_size;
            noninert_amean(:,j+1) = noninert_amean(:,j+1)+noninert_ensemble(4*k-3:4*k,j+1)/ens_size;
            order1_amean(:,j+1) = order1_amean(:,j+1)+order1_ensemble(4*k-3:4*k,j+1)/ens_size;
        end
        
        avg_perf_error(j+1) = avg_perf_error(j+1)+(1/num_trials)*...
            sqrt(((perf_amean(1,j+1)-truth(1,time_rat*j+1))^2+(perf_amean(2,j+1)-truth(2,time_rat*j+1))^2));
        avg_noninert_error(j+1) = avg_noninert_error(j+1)+(1/num_trials)*...
            sqrt(((noninert_amean(1,j+1)-truth(1,time_rat*j+1))^2+(noninert_amean(2,j+1)-truth(2,time_rat*j+1))^2));
        avg_order1_error(j+1) = avg_order1_error(j+1)+(1/num_trials)*...
            sqrt(((order1_amean(1,j+1)-truth(1,time_rat*j+1))^2+(order1_amean(2,j+1)-truth(2,time_rat*j+1))^2));
        
    end
    
end

figure
plot(avg_obs_error)
hold on
plot(avg_perf_error,'r')
plot(avg_noninert_error,'k')
plot(avg_order1_error,'g')
legend('Observational Error','Perfect Model Error','Noninertial Model Error','Order 1 Model Error')
title('Average RMS Error for Different Models')
xlabel('Time')
ylabel('Avg RMS Error')