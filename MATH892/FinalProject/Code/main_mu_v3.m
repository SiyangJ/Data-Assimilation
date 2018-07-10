mus = logspace(-4,0.6,10);
nmu = length(mus);
Afacs = [0,0.2];
Bfacs = [0,0.2];
nens = 40;
deltaobs = 0.05;

%%
BMss = cell(1,4);
for j=1:4
    
    if j==1
        A = 0.2; B = 0.2;
    elseif j==2
        A = 0; B = 0.2;
    elseif j==3
        A = 0.2; B = 0;
    elseif j==4
        A = 0; B = 0;
    end
    
    
BMss{j} = cell(nmu,1);

for i=1:nmu
    BMss{j}{i} = EnKFfun_v3(mus(i),A,B,nens,deltaobs);
end

fprintf('Completed BMss{%d}\n',j)
end
save('data/mutest2.mat','BMss')

%%

figure
for j = 1:4

BMs = BMss{j};

eaps = zeros(1,nmu);
eab1s = zeros(1,nmu);
eab2s = zeros(1,nmu);
eab3s = zeros(1,nmu);

for i=1:nmu
    
    eaps(i) = BMs{i}{6};
    eab1s(i) = BMs{i}{7};
    eab2s(i) = BMs{i}{8};
    eab3s(i) = BMs{i}{9};

end

subplot(2,2,j)
hold on
plot(mus,eaps,'-o')
plot(mus,eab1s,'-s')
plot(mus,eab2s,'-*')
plot(mus,eab3s,'-d')
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend({'No Bias Model','Bias Model I','Bias Model II','Bias Model III'})
xlabel('\mu')
ylabel('$$\langle\langle\bar{e_{a}}\rangle\rangle$$','interpreter','latex')
title(sprintf('A=%f, B=%f',BMs{1}{2},BMs{1}{3}))
end
suptitle(sprintf('deltaobs=%f, nens=%f',BMss{1}{1}{4},BMss{1}{1}{5}))
