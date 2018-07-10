%%
figure
for m = 0:2
for i=1:4
    subplot(3,4,m*4+i)
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
        ami = xamb2(i,:);
        fmi = xfmb2(i,:);
        ame = reshape(xacovb2(i,i,:),[1,nobs+1]);
        fme = reshape(xfcovb2(i,i,:),[1,nobs+1]);
    end
    errorbar(tobs,ami,ame,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    errorbar(tobs,fmi,fme,'r-s','LineWidth',2,'MarkerSize',8)
    plot(tobs,xmf(i,:),'g--*')
    %xlim([10,15])
    hold off
    title(sprintf('Model %d, Variable %d',m,i))
end
end

%%

tspan = floor(nobs/8)+1:nobs+1;

errpmat = xfm(1:ndim,tspan)-truestate(1:ndim,tspan);
errp = mean(errpmat,2);
errpcov = diag(cov(errpmat')).^0.5;
errb1mat = xfmb1(1:ndim,tspan)-truestate(1:ndim,tspan);
errb1 = mean(errb1mat,2);
errb1cov = diag(cov(errb1mat')).^0.5;
errb2mat = xfmb2(1:ndim,tspan)-truestate(1:ndim,tspan);
errb2 = mean(errb2mat,2);
errb2cov = diag(cov(errb2mat')).^0.5;

%%
figure
hold on

errorbar(errp,errpcov)
errorbar(errb1,errb1cov)
errorbar(errb2,errb2cov)

hold off

legend({'perfect','linear','order1'})

%%
sumrmsap = 0;
sumrmsab1 = 0;
sumrmsab2 = 0;

tspan = floor(nobs/2)+1:nobs+1;

%tspan = 50:nobs+1;
t1 = tspan(1);
tf = tspan(end);
t0 = tspan(1)-1;
T = length(tspan);

% Calculat the change in the time averaged analysis error
% To see the convergence of the time averaged error
% Useful for inspection of settling time

eaparr = zeros(T,1);
eab1arr = zeros(T,1);
eab2arr = zeros(T,1);

for t = tspan
    
    sumrmsap = sumrmsap + norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    eaparr(t-t0) = sumrmsap/(t-t0);
    sumrmsab1 = sumrmsab1 + norm(xamb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    eab1arr(t-t0) = sumrmsab1/(t-t0);
    sumrmsab2 = sumrmsab2 + norm(xamb2(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    eab2arr(t-t0) = sumrmsab2/(t-t0);
    
end

figure
hold on
fplot(@(x)eap,[t1,tf],'--r','LineWidth',2)
fplot(@(x)eab1,[t1,tf],'--g','LineWidth',2)
fplot(@(x)eab2,[t1,tf],'--b','LineWidth',2)

plot(tspan,eaparr,'-r')
plot(tspan,eab1arr,'-g')
plot(tspan,eab2arr,'-b')

hold off
legend({'AvP','AvB1','AvB2','p','b1','b2'})

%%
%tspan = floor(nobs/2)+1:nobs+1;

figure
for m = 1:10
    
load(sprintf('data/LSWLongRunMu%d.mat',m))

tspan = 1:nobs+1;
t1 = tspan(1);
tf = tspan(end);
t0 = tspan(1)-1;
T = length(tspan);

% Calculat the change in the time averaged analysis error
% To see the convergence of the time averaged error
% Useful for inspection of settling time

eapts = zeros(T,1);
eab1ts = zeros(T,1);
eab2ts = zeros(T,1);

for t = tspan
    
    eapts(t-t0) = norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    eab1ts(t-t0) = norm(xamb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    eab2ts(t-t0) = norm(xamb2(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    
end

subplot(5,2,m)
hold on
fplot(@(x)eap,[t1,tf],'--r','LineWidth',2)
fplot(@(x)eab1,[t1,tf],'--g','LineWidth',2)
fplot(@(x)eab2,[t1,tf],'--b','LineWidth',2)

plot(tspan,eapts,'-r')
plot(tspan,eab1ts,'-g')
plot(tspan,eab2ts,'-b')

hold off
%legend({'AvP','AvB1','AvB2','p','b1','b2'})
title(sprintf('RMS at time frame, mu = %f',inflmu))
set(gca,'YScale','log')

end
