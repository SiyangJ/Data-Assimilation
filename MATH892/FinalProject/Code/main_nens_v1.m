nenss = 40:10:100;
nnens = length(nenss);
% mus chosen for best in each situation
% Need to investigate further Type C model error,
% which is inconsistent with the paper
% 4 is best for all models
% Doesn't make much sense
mus = [1.23,0.2782,0.3782,0.3782];
Afacs = [0,0.2];
Bfacs = [0,0.2];

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
    
    
BMss{j} = cell(nnens,1);

for i=1:nnens
    BMss{j}{i} = EnKFfun_v2(mus(j),A,B,nenss(i));
end

%save('data/mutest1.mat','BMs_mu_2','-append')
fprintf('Completed BMss{%d}\n',j)
end
save('data/nenstest1.mat','BMss')

%%

figure
for j = 1:4

BMs = BMss{j};

eaps = zeros(1,nnens);
eab1s = zeros(1,nnens);
eab2s = zeros(1,nnens);
eab3s = zeros(1,nnens);

for i=1:nnens
    
    eaps(i) = BMs{i}{6};
    eab1s(i) = BMs{i}{7};
    eab2s(i) = BMs{i}{8};
    eab3s(i) = BMs{i}{9};

end

subplot(2,2,j)
hold on
plot(nenss,eaps,'-o')
plot(nenss,eab1s,'-s')
plot(nenss,eab2s,'-*')
plot(nenss,eab3s,'-d')
hold off
%set(gca,'XScale','log')
%set(gca,'YScale','log')
legend({'No Bias Model','Bias Model I','Bias Model II','Bias Model III'})
xlabel('Ensemble Size')
ylabel('$$\langle\langle\bar{e_{a}}\rangle\rangle$$','interpreter','latex')
title(sprintf('A=%f, B=%f',BMs{1}{2},BMs{1}{3}))
end
suptitle(sprintf('deltaobs=%f, mu=%f',BMss{1}{1}{4},BMss{1}{1}{1}))
