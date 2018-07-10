mus = logspace(-4,0.6,10);
nmu = length(mus);

epsilon = 0.01;

%%
BMs = cell(1,nmu);

for i=1:nmu
    BMs{i} = EnKFfun_LSW_v1(mus(i),epsilon,sprintf('data/LSWLongRunMu%d.mat',i));
end

save('data/LSWmutest1.mat','BMs')

%%

figure

eaps = zeros(1,nmu);
eab1s = zeros(1,nmu);
eab2s = zeros(1,nmu);

for i=1:nmu
    
    eaps(i) = BMs{i}{5};
    eab1s(i) = BMs{i}{6};
    eab2s(i) = BMs{i}{7};

end

hold on
plot(mus,eaps,'-o')
plot(mus,eab1s,'-s')
plot(mus,eab2s,'-*')
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend({'Perfect Model','Linear Model','Order 1 Inertia'})
xlabel('\mu')
ylabel('$$\langle\langle\bar{e_{a}}\rangle\rangle$$','interpreter','latex')
