mus = logspace(-4,0.6,10);
nmu = length(mus);

gamma = 0.05;

%%
BMs = cell(1,nmu);

for i=1:nmu
    BMs{i} = EnKFfun_ME2_v1(mus(i),gamma,sprintf('data/ME2LongRunMu%d.mat',i));
end

save('data/ME2mutest1.mat','BMs')

%%

figure

eaps = zeros(1,nmu);
eab1s = zeros(1,nmu);
%eab2s = zeros(1,nmu);
eab3s = zeros(1,nmu);

for i=1:nmu
    
    eaps(i) = BMs{i}{7};
    eab1s(i) = BMs{i}{8};
    %eab2s(i) = BMs{i}{8};
    eab3s(i) = BMs{i}{9};

end

hold on
plot(mus,eaps,'-o')
plot(mus,eab1s,'-s')
%plot(mus,eab2s,'-*')
plot(mus,eab3s,'-d')
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend({'No Bias Model','Bias Model I','Bias Model III'})
xlabel('\mu')
ylabel('$$\langle\langle\bar{e_{a}}\rangle\rangle$$','interpreter','latex')
