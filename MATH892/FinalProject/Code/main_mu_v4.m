mus = logspace(-4,0.6,10);
nmu = length(mus);
Afacs = [0,0.2];
Bfacs = [0,0.2];

%%
% BMs_mu_1 = cell(nmu,1);
% 
% for i=1:nmu
%     BMs_mu_1{i} = EnKFfun_v1(mus(i),0.2,0.2);
% end
% 
% save('data/mutest1.mat','BMs_mu_1')
%%

BMs_mu_2 = cell(nmu,1);

for i=1:nmu
    BMs_mu_2{i} = EnKFfun_v1(mus(i),0,0.2);
end

save('data/mutest1.mat','BMs_mu_2','-append')
fprintf('Saved BMs_mu_2\n')
%%

BMs_mu_3 = cell(nmu,1);

for i=1:nmu
    BMs_mu_3{i} = EnKFfun_v1(mus(i),0.2,0);
end

save('data/mutest1.mat','BMs_mu_3','-append')
fprintf('Saved BMs_mu_3\n')

%%
BMs_mu_4 = cell(nmu,1);

for i=1:nmu
    BMs_mu_4{i} = EnKFfun_v1(mus(i),0,0);
end

save('data/mutest1.mat','BMs_mu_4','-append')
fprintf('Saved BMs_mu_4\n')

%%

BMss = {BMs_mu_1,BMs_mu_2,BMs_mu_3,BMs_mu_4};

%figure
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

%subplot(2,2,j)
figure
hold on
plot(mus,eaps,'-o')
plot(mus,eab1s,'-s')
plot(mus,eab2s,'-*')
plot(mus,eab3s,'-d')
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend({'No Bias Estimation','Bias Model I','Bias Model II','Bias Model III'})
xlabel('\mu','FontSize',20)
ylabel('$$\langle\langle\bar{\mathbf{e}}_a\rangle\rangle$$','interpreter','latex','FontSize',20)
title(sprintf('A=%f, B=%f',BMs{1}{2},BMs{1}{3}),'FontSize',16,'FontWeight','normal')
end
%suptitle(sprintf('deltaobs=%f, nens=%f',BMss{1}{1}{4},BMss{1}{1}{5}))
