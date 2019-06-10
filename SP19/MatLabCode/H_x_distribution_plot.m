clc;clear;close all

xrange = [-15,0,15];
yrange = [-15,0,15];

for ix = 1:length(xrange)
    for iy = 1:length(yrange)
%%
    subplot(length(xrange),length(yrange),(ix-1)*length(yrange) + iy)
    nens = 10e6;
    ndim = 2;

    V = [1,-1;1,1];
    S = [10,0;0,1];

    M = V'*S*V;

    ens = mvnrnd([xrange(ix),yrange(iy)],[3.2,0;0,3.2],nens);
    ens_mean = mean(ens,1);
    fprintf('X mean is:')
    disp(ens_mean)

    % ens(ens(:,2)<0) = 0;
    %%
%     figure
%     nbins=[50 50];
%     [N,C]=hist3(ens,nbins);
%     lN = log(N);
%     lN(lN==-inf)=0;
%     contourf(C{1},C{2},lN',50)

    %%

    H = [0,1;1,0];

    %h = @(x) EW09_obs(x')';
    %h = @(x) [x(:,1) + x(:,2),x(:,1)];
    h = @(x) [(x(:,1)+x(:,2)).^2,(x(:,1)+x(:,2)).^3];
    
%     yensp = h(ens);
%     yens  = yensp(:,[1,4]);

    yens = h(ens);
    
    h_ens_mean = h(ens_mean);
    fprintf('h(X mean) is:')
    disp(h_ens_mean)
    
    yens_mean = mean(yens,1);
    fprintf('Y mean is:')
    disp(yens_mean)

    %%
    nbins=[10 10];
    [N,C]=hist3(yens,nbins);
    lN = log(N);
    lN(lN==-inf)=0;
    %%
%     figure
    contourf(C{1},C{2},lN',10)
    hold on
    scatter(h_ens_mean(:,1),h_ens_mean(:,2),'*')
    scatter(yens_mean(:,1),yens_mean(:,2),'o')
    end
end
