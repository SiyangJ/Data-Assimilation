clear all;
syms  x1 x2;   %variables must be named x1...xn
	f=[x1+.1*x2+.005;x2+.1];
	h=[x1];
	x_tr=[1;1]; %initial value of state
    num_members=10; % number of ensemble members to take
    w=[10^-3; .02];  %process noise standard deviation
    z=[10];  %measurement noise standard deviation
    ICN=[10; 10]; % inital condition standard deaviation, error in knowing where to start
    x_ini=ones(2,num_members)+ICN.*randn(2,num_members);
	num_iterations=20;

[truth,obs,estimate]=ensemblekfilter(f,h,x_tr,x_ini,w,z,num_iterations);
    
plot(truth(1,:))
hold on
display('This is the true dynamics for x1 Hit enter to continue plotting');
pause;
plot(obs(1,:))
display('These are the observations of x1. Hit enter to continue plotting' );
pause;
hold on
plot(estimate(1,:))
display(' These are the analysis estimates for x1. Hit enter to continue plotting');
pause;

figure

plot(truth(2,:))
display('This is the true dynamics for x2 Hit enter to continue plotting');
hold on
plot(estimate(2,:))
display('These are the analysis estimates x2');


