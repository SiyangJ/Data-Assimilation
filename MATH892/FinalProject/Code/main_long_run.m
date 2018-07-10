%%
% Test settling time for long run
% 800 seems to be enough
% Test for 1000
% initial spread enlarged


%%
deltaobs = 0.05;
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

%%
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
    
    EnKFfun_v4(mus(j),A,B,nens,deltaobs,sprintf('data/AllDataLongRunCase%d.mat',j));

    fprintf('Completed BMss{%d}\n',j)
end

%%




