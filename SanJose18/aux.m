%%

for m = [0,1,3]
figure
for i=1:4
    subplot(1,4,i)
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
end
suptitle(sprintf('Model %d',m))
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

tspan = 101:201;

zetab1 = mean(xfmb1(ndim+1:end,tspan),2);
zetab1cov = diag(cov(xfmb1(ndim+1:end,tspan)')).^0.5;
xib2 = mean(xfmb2(ndim+1:end,tspan),2);
xib2cov = diag(cov(xfmb2(ndim+1:end,tspan)')).^0.5;
zetab3 = mean(xfmb3(ndim+1:2*ndim,tspan),2);
zetab3cov = diag(cov(xfmb3(ndim+1:2*ndim,tspan)')).^0.5;
xib3 = mean(xfmb3(2*ndim+1:end,tspan),2);
xib3cov = diag(cov(xfmb3(2*ndim+1:end,tspan)')).^0.5;
%%
%xizeta = 0;
figure
subplot(2,1,1)
hold on
plot(zeta*deltaobs,'b-','LineWidth',2)
plot(zetab1,'k--*')
plot(zetab3,'r--o','MarkerSize',7,'MarkerEdgeColor','r')
%errorbar(zetab1,zetab1cov)
%errorbar(zetab3,zetab3cov)
%errorbar(xib2,xib2cov)
hold off
legend({'Truth','Bias Model 1','Bias Model 3'},'FontSize',12)
xlabel('i','FontSize',16)
ylabel('$$\langle\bar{\mathbf{b}}_{a}^{(i)}\rangle$$','interpreter','latex',...
   'FontSize',20)
subplot(2,1,2)
hold on
plot(-xi,'b-','LineWidth',2)
plot(xib2,'k--*')
plot(xib3,'r--o','MarkerSize',7,'MarkerEdgeColor','r')
hold off
legend({'Truth','Bias Model 2','Bias Model 3'},'FontSize',12,'Location','northwest')   
xlabel('i','FontSize',16)
ylabel('$$\langle\bar{\mathbf{c}}_{a}^{(i)}\rangle$$','interpreter','latex',...
    'FontSize',20)
%suptitle('Estimation of Bias, Type III Model Error','FontSize',16)

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
sumrmsp = 0;
sumrmsb1 = 0;
sumrmsb2 = 0;
sumrmsb3 = 0;

tspan = 51:100;
T = length(tspan);

ep = sumrmsp/T;
eb1 = sumrmsb1/T;
eb2 = sumrmsb2/T;
eb3 = sumrmsb3/T;

for t = tspan
    
    sumrmsp = sumrmsp + norm(xfm(:,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb1 = sumrmsb1 + norm(xfmb1(1:ndim,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb2 = sumrmsb2 + norm(xfmb2(1:ndim,t)+xfmb2(ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    sumrmsb3 = sumrmsb3 + norm(xfmb3(1:ndim,t)+xfmb3(2*ndim+1:end,t)-truestate(:,t))/sqrt(ndim);
    
end

ep = sumrmsp/T;
eb1 = sumrmsb1/T;
eb2 = sumrmsb2/T;
eb3 = sumrmsb3/T;

