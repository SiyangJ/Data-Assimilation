deltaobss = 0.03:0.01:0.12;
ndeltaobs = length(deltaobss);
% mus chosen for best in each situation
% Need to investigate further Type C model error,
% which is inconsistent with the paper
% 4 is best for all models
% Doesn't make much sense
mus = [1.23,0.3782,0.3782,0.3782];
Afacs = [0,0.2];
Bfacs = [0,0.2];
% nens = 60 seems to be a good choice.
nens = 60;

xlog = false;
ylog = true;

%%
BMss = cell(1,4);
for j=1:4
    
    if j==1
        A = 0.2; B = 0.2;
    elseif j==2
        A = 0;   B = 0.2;
    elseif j==3
        A = 0.2; B = 0;
    elseif j==4
        A = 0;   B = 0;
    end
    
    
BMss{j} = cell(ndeltaobs,1);

for i=1:ndeltaobs
    BMss{j}{i} = EnKFfun_v3(mus(j),A,B,nens,deltaobss(i));
end

fprintf('Completed BMss{%d}\n',j)
end
save('data/deltaobstest1.mat','BMss')

%%

figure
for j = 1:4

BMs = BMss{j};

eaps = zeros(1,ndeltaobs);
eab1s = zeros(1,ndeltaobs);
eab2s = zeros(1,ndeltaobs);
eab3s = zeros(1,ndeltaobs);

for i=1:ndeltaobs
    
    eaps(i) = BMs{i}{6};
    eab1s(i) = BMs{i}{7};
    eab2s(i) = BMs{i}{8};
    eab3s(i) = BMs{i}{9};

end

subplot(2,2,j)
hold on
plot(deltaobss,eaps,'-o')
plot(deltaobss,eab1s,'-s')
plot(deltaobss,eab2s,'-*')
plot(deltaobss,eab3s,'-d')
hold off
if xlog
    set(gca,'XScale','log')
end
if ylog
    set(gca,'YScale','log')
end
legend({'No Bias Model','Bias Model I','Bias Model II','Bias Model III'})
xlabel('Delta Observations')
ylabel('$$\langle\langle\bar{e_{a}}\rangle\rangle$$','interpreter','latex')
title(sprintf('A=%f, B=%f',BMs{1}{2},BMs{1}{3}))
end
suptitle(sprintf('Analysis Error, nens=%d',nens))

%%

figure
for j = 1:4

BMs = BMss{j};

eps = zeros(1,ndeltaobs);
eb1s = zeros(1,ndeltaobs);
eb2s = zeros(1,ndeltaobs);
eb3s = zeros(1,ndeltaobs);

for i=1:ndeltaobs
    
    eps(i) = BMs{i}{10};
    eb1s(i) = BMs{i}{11};
    eb2s(i) = BMs{i}{12};
    eb3s(i) = BMs{i}{13};

end

subplot(2,2,j)
hold on
plot(deltaobss,eps,'-o')
plot(deltaobss,eb1s,'-s')
plot(deltaobss,eb2s,'-*')
plot(deltaobss,eb3s,'-d')
hold off
if xlog
    set(gca,'XScale','log')
end
if ylog
    set(gca,'YScale','log')
end
legend({'No Bias Model','Bias Model I','Bias Model II','Bias Model III'})
xlabel('Delta Observations')
ylabel('$$\langle\langle\bar{e_{f}}\rangle\rangle$$','interpreter','latex')
title(sprintf('A=%f, B=%f',BMs{1}{2},BMs{1}{3}))
end
suptitle(sprintf('Forecast Error, nens=%d',nens))
