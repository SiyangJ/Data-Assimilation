clear;set(0,'defaultaxesfontsize',20);format long
%%% p7.m weak 4DVAR for sin map (Ex. 1.3)
%% setup

%Vary these parameters:
J=10;% number of steps
gamma=1e0;% observational noise variance is gamma^2
C0=1;% prior initial condition variance
m0=0;% prior initial condition mean
sd=1;rng(sd);% choose random number seed

%% truth
alpha=2.5;% dynamics determined by alpha

vt(1)=sqrt(C0)*randn;% truth initial condition
for j=1:J
    vt(j+1)=alpha*sin(vt(j)); % create truth
    y(j)=vt(j+1)+gamma*randn;% create observations 
end

%% Generate initial guess and solve optimization problem 

    uu=-0.1;% initial guess  
   %uu=3;

% solve with matlab fminsearch with defaults 
% exitflag=1 ==> convergence
%%%%%%
op = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','dfp');
op = optimoptions(op,'StepTolerance',1e-15,'OptimalityTolerance',1e-11);
[uin,fval,exitflag]=fminunc(@(u)CostFunction(u,y,gamma,alpha,m0,C0,J),uu,op)
%%%%%%
%op = optimset('TolFun',1e-10);
%[uin,fval,exitflag]=fminsearch(@(u)CostFunction(u,y,gamma,alpha,m0,C0,J),uu,op)
%%%%%%

%Generate map from solution uin and from initial guess uu 
vmap(1)=uin;
trajig(1)=uu;
for j=2:J+1
vmap(j)=alpha*sin(vmap(j-1));
trajig(j)=alpha*sin(trajig(j-1));
end

%Plot
figure;plot([0:J],vmap,'Linewidth',2);hold;plot([0:J],vt,'r','Linewidth',2)
plot([0:J],trajig,'k','Linewidth',2);
plot([1:J],y,'g','Linewidth',2);
hold;xlabel('j');legend('Map','Truth','Init Guess','Obs','Location','northwest')