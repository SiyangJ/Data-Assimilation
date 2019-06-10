clear; set(0,'defaultaxesfontsize',20); format long
%%% p3.m MCMC RWM algorithm for logistic map (Ex. 1.4)
%% setup
J=5;% number of steps
r=4;% dynamics determined by alpha
gamma=0.2;% observational noise variance is gamma?2
C0=0.01;% prior initial condition variance
m0=0.5;% prior initial condition mean
sd=10;rng(sd);% choose random number seed

%% truth

vt=0.3;vv(1)=vt;% truth initial condition
Jdet=1/2/C0*(vt-m0)^2;% background penalization
Phidet=0;% initialization model-data misfit functional
for j=1:J
% can be replaced by Psi for each problem
vv(j+1)=r*vv(j)*(1-vv(j));% create truth
y(j)=vv(j+1)+gamma*randn;% create data
Phidet=Phidet+1/2/gamma^2*(y(j)-vv(j+1))^2;% misfit functional
end
Idet=Jdet+Phidet;% compute log posterior of the truth

%% solution
% Markov Chain Monte Carlo: N forward steps of the
% Markov Chain on R (with truth initial condition)
N=1e5;% number of samples
V=zeros(N,1);% preallocate space to save time
beta=0.05;% step-size of random walker
v=vt;% truth initial condition (or else update I0)
n=1; bb=0; rat(1)=0;
while n<=N
w=v+sqrt(2*beta)*randn;% propose sample from random walker
vv(1)=w;
Jdetprop=1/2/C0*(w-m0)^2;% background penalization
Phidetprop=0;
for i=1:J
vv(i+1)=r*vv(i)*(1-vv(i));
Phidetprop=Phidetprop+1/2/gamma^2*(y(i)-vv(i+1))^2;
end
Idetprop=Jdetprop+Phidetprop;% compute log posterior of the proposal

if rand<exp(Idet-Idetprop)% accept or reject proposed sample
v=w; Idet=Idetprop; bb=bb+1;% update the Markov chain
end
rat(n)=bb/n;% running rate of acceptance
V(n)=v;% store the chain
n=n+1;
end
dx=0.0005; v0=[0.01:dx:0.99];
Z=hist(V,v0);% construct the posterior histogram
figure(1), plot(v0,Z/trapz(v0,Z),'k','Linewidth',2)% visualize the posterior