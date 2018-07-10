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
