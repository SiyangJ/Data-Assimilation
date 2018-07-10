
%%
figure
for m = 0:3
%for i=1:4
    subplot(4,1,m+1)
    i = 10;
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
        ami = xamb2(i,:)+xamb2(ndim+i,:);
        fmi = xfmb2(i,:)+xfmb2(ndim+i,:);
        ame = reshape(xacovb2(i,i,:),[1,nobs+1]);
        fme = reshape(xfcovb2(i,i,:),[1,nobs+1]);
    elseif m==3
        ami = xamb3(i,:)+xamb3(2*ndim+i,:);
        fmi = xfmb3(i,:)+xfmb3(2*ndim+i,:);
        ame = reshape(xacovb3(i,i,:),[1,nobs+1]);
        fme = reshape(xfcovb3(i,i,:),[1,nobs+1]);
    end
    errorbar(tobs,ami,ame,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    errorbar(tobs,fmi,fme,'r-s','LineWidth',2,'MarkerSize',8)
    plot(tobs,xmf(i,:),'g--*')
    hold off
%end
end

%%
for m = 1:2
figure
for i=1:4
    subplot(1,4,i)
    hold on
    %plot(ttraj,truetraj(i,:),'b-','LineWidth',4)
    fplot(@(x)zeta(i)*deltaobs,[0,tend],'b-','LineWidth',4)
    if m==1
        ami = xamb1(ndim+i,:);
        fmi = xfmb1(ndim+i,:);
        ame = reshape(xacovb1(ndim+i,ndim+i,:),[1,nobs+1]).^0.5;
        fme = reshape(xfcovb1(ndim+i,ndim+i,:),[1,nobs+1]).^0.5;
    elseif m==2
        ami = xamb3(ndim+i,:);
        fmi = xfmb3(ndim+i,:);
        ame = reshape(xacovb3(ndim+i,ndim+i,:),[1,nobs+1]).^0.5;
        fme = reshape(xfcovb3(ndim+i,ndim+i,:),[1,nobs+1]).^0.5;
    end
    errorbar(tobs,ami,ame,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    errorbar(tobs,fmi,fme,'r-s','LineWidth',2,'MarkerSize',8)
    %plot(tobs,xmf(i,:),'g--*')
    hold off
end
end

%%
figure
hold on
errorbar(zetab1,zetab1cov)
%errorbar(xib2,xib2cov)
plot(zeta*deltaobs)
%plot(-xi)
hold off

%%

tspan = floor(nobs/8)+1:nobs+1;

errpmat = xfm(1:ndim,tspan)-truestate(1:ndim,tspan);
errp = mean(errpmat,2);
errpcov = diag(cov(errpmat')).^0.5;
errb1mat = xfmb1(1:ndim,tspan)-truestate(1:ndim,tspan);
errb1 = mean(errb1mat,2);
errb1cov = diag(cov(errb1mat')).^0.5;
errb2mat = xfmb2(1:ndim,tspan)+xfmb2(ndim+1:end,tspan)-truestate(1:ndim,tspan);
errb2 = mean(errb2mat,2);
errb2cov = diag(cov(errb2mat')).^0.5;
errb3mat = xfmb3(1:ndim,tspan)+xfmb3(2*ndim+1:end,tspan)-truestate(1:ndim,tspan);
errb3 = mean(errb3mat,2);
errb3cov = diag(cov(errb3mat')).^0.5;

%%
figure
hold on

errorbar(errp,errpcov)
errorbar(errb1,errb1cov)
errorbar(errb2,errb2cov)
errorbar(errb3,errb3cov)

hold off

legend({'p','b1','b2','b3'})

%%
sumrmsap = 0;
sumrmsab1 = 0;
sumrmsab2 = 0;
sumrmsab3 = 0;

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
eab3arr = zeros(T,1);

for t = tspan
    
    sumrmsap = sumrmsap + norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    eaparr(t-t0) = sumrmsap/(t-t0);
    sumrmsab1 = sumrmsab1 + norm(xamb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    eab1arr(t-t0) = sumrmsab1/(t-t0);
    sumrmsab2 = sumrmsab2 + norm(xamb2(1:ndim,t)+xamb2(ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    eab2arr(t-t0) = sumrmsab2/(t-t0);
    sumrmsab3 = sumrmsab3 + norm(xamb3(1:ndim,t)+xamb3(2*ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    eab3arr(t-t0) = sumrmsab3/(t-t0);
    
end

figure
hold on
fplot(@(x)eap,[t1,tf],'--r','LineWidth',2)
fplot(@(x)eab1,[t1,tf],'--g','LineWidth',2)
fplot(@(x)eab2,[t1,tf],'--b','LineWidth',2)
fplot(@(x)eab3,[t1,tf],'--c','LineWidth',2)

plot(tspan,eaparr,'-r')
plot(tspan,eab1arr,'-g')
plot(tspan,eab2arr,'-b')
plot(tspan,eab3arr,'-c')

hold off
legend({'AvP','AvB1','AvB2','AvB3','p','b1','b2','b3'})

%%
%tspan = floor(nobs/2)+1:nobs+1;

for m = 1:10
    
load(sprintf('data/ME2LongRunMu%d.mat',m))

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
%eab2ts = zeros(T,1);
eab3ts = zeros(T,1);

for t = tspan
    
    eapts(t-t0) = norm(xam(:,t)-truestate(:,t))/sqrt(ndim);
    eab1ts(t-t0) = norm(xamb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    %eab2ts(t-t0) = norm(xamb2(1:ndim,t)+xamb2(ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    eab3ts(t-t0) = norm(xamb3(1:ndim,t)+xamb3(2*ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    
end

figure
hold on
fplot(@(x)eap,[t1,tf],'--r','LineWidth',2)
fplot(@(x)eab1,[t1,tf],'--g','LineWidth',2)
%fplot(@(x)eab2,[t1,tf],'--b','LineWidth',2)
fplot(@(x)eab3,[t1,tf],'--c','LineWidth',2)

plot(tspan,eapts,'-r')
plot(tspan,eab1ts,'-g')
%plot(tspan,eab2ts,'-b')
plot(tspan,eab3ts,'-c')

hold off
legend({'AvP','AvB1','AvB3','p','b1','b3'})
title(sprintf('RMS at time frame, mu = %f',inflmu))
set(gca,'YScale','log')

end
