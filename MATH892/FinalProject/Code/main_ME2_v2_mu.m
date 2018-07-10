mus = [1e-4,3e-4,7e-4,1e-3,3e-3,7e-3,1e-2,5e-2,1e-1,5e-1,1,2,4,6];
nmu = length(mus);

gamma = 0.05;

%%
BMs = cell(1,nmu);

for i=1:nmu
    BMs{i} = EnKFfun_ME2_v2(mus(i),gamma,sprintf('data/ME2v2LongRunMu%d.mat',i));
end

save('data/ME2v2mutest1.mat','BMs')

%%

figure

eaps = zeros(1,nmu);
eab1s = zeros(1,nmu);

for i=1:nmu
    
    eaps(i) = BMs{i}{5};
    eab1s(i) = BMs{i}{6};

end

hold on
plot(mus,eaps,'-o')
plot(mus,eab1s,'-s')
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend({'No Bias Model','Bias Model'})
xlabel('\mu')
ylabel('$$\langle\langle\bar{e_{a}}\rangle\rangle$$','interpreter','latex')
