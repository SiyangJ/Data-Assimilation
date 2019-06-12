%% A simple 2-dimensional case

for ix = 1:2
    
for iy = 1:2

ndim = 2;
dobs = 2;

Pf = 1.5*eye(ndim);
truestate = [4,2]';
xinit = mvnrnd(truestate,Pf)';

sigmaobs = 0.16;
Robs = sigmaobs.^2.* eye(dobs);

% H = [1,1];
H = eye(2);

yobs = mvnrnd(H * truestate, Robs,1)';
disp(yobs)
s = sum(yobs);
yobs(1) = yobs(1) * (0.5 + 0.01*randn);
yobs(2) = s - yobs(1);
disp(yobs)

K = Pf * H' / (H * Pf * H' + Robs);

xpost = xinit + K * (yobs - H * xinit);
Pa = (eye(2) - K * H) * Pf;
%Pa = (Pa + Pa')/2;

nens = 1e7;

ens = mvnrnd(xpost,Pa,nens);

subplot(2,2,(ix-1)*2+iy)
nbins=[100 100];
[N,C]=hist3(ens,nbins);
lN = log(N);
lN(lN==-inf)=0;
contourf(C{1},C{2},lN',8)
% contourf(C{1},C{2},N',10)
hold on
scatter(xinit(1),xinit(2),'o')
scatter(truestate(1),truestate(2),'*')
scatter(xpost(1),xpost(2),'p')
hold off
end
end

%% A 2-dimensional case with particle filter

PLOT = false;

inflmuarray = [0.05,0.1,0.2,0.4,0.8];
n_ensarray = [10,20,100,1e4];
ixarray = 1:100;

errorarray = zeros(length(inflmuarray),length(n_ensarray),length(ixarray),5);

for i_inflmu = 1:length(inflmuarray)
inflmu = inflmuarray(i_inflmu);
for i_n_ens = 1:length(n_ensarray)
n_ens = n_ensarray(i_n_ens);
if PLOT
figure('NumberTitle', 'off', 'Name', sprintf('inflmu=%3.2f, nens=%d',inflmu,n_ens));
end
for ix = ixarray
ndim = 2;
dobs = 2;

% Truth
truestate = xlong(mod(2*ix,30)+10:mod(2*ix,30)+11,mod(127*ix,10000)+1);

% Obs
sigmaobs = 0.16;
Robs = sigmaobs.^2.* eye(dobs);

% H = [1,1];
H = eye(2);

h = @(x) H * x;%@two_pair_stupid;

% hidden_obs = H * truestate;
yobs = mvnrnd(h(truestate), Robs,1)';

% Initial
Pf = 1.5*eye(ndim);
xinit = mvnrnd(truestate,Pf)';

% Ensemble/Particle
xens = mvnrnd(xinit,Pf,n_ens)';
pens = mvnrnd(xinit,Pf,n_ens)';
W = ones(n_ens,1)/n_ens;

pf = cov(xens');
pf = pf+inflmu*trace(pf)/ndim*eye(ndim);
pfht = pf*H';   
K = pfht*pinv(H*pfht+Robs);  

xenspost = xens + K * (mvnrnd(yobs,Robs,n_ens)' - H * xens);

xpost = mean(xenspost,2);
Pa = cov(xenspost');

innov = yobs - h(pens);

% Jack's original code
% Wtmp =  -0.5*innov'/Robs*innov;    
% Wmax=max(Wtmp);
% Wtmp  = Wtmp-Wmax;
% W = W.*exp(Wtmp'); 
% W=W/sum(W);

% My correction to his code
Wtmp =  -0.5*diag(innov'/Robs*innov);    
Wmax=max(Wtmp);
Wtmp  = Wtmp-Wmax;
W = W.*exp(Wtmp); 
W=W/sum(W);

resamp_thresh = 0.5;
wiggle = 0.1;

ppost1 = pens * W;

if   1/sum(W.^2)/n_ens < resamp_thresh %test for resampling
    fprintf('Resampling triggered, effective ratio is %6.4f.\n',1/sum(W.^2)/n_ens)
    sampIndex = resampleMultinomial(W);
    pens = pens(:,sampIndex)+wiggle*randn(ndim,n_ens);
    W = ones(n_ens,1)/n_ens;
end

ppost = pens * W;

nens = 1e7;

if PLOT

% Partical Filter
subplot(2,6,ix)
nbin = min(floor(n_ens/5),50);
nbins=[nbin,nbin];
[N,C]=hist3w(pens',W,nbins);
lN = log(N);
lN(lN==-inf)=0;
nlevel = min(nbin,6);
contourf(C{1},C{2},lN',nlevel)
% contourf(C{1},C{2},N',10)
hold on
scatter(xinit(1),xinit(2),'o')
scatter(truestate(1),truestate(2),'*')
scatter(ppost(1),ppost(2),'p')
hold off
title(sprintf('truestate=[%4.3f,%4.3f]',truestate(1),truestate(2)))

% EnKF
subplot(2,6,ix + 6)
nbin = min(floor(n_ens/5),50);
nbins=[nbin,nbin];
[N,C]=hist3(xenspost',nbins);
lN = log(N);
lN(lN==-inf)=0;
nlevel = min(nbin,6);
contourf(C{1},C{2},lN',nlevel)
% contourf(C{1},C{2},N',10)
hold on
scatter(xinit(1),xinit(2),'o')
scatter(truestate(1),truestate(2),'*')
scatter(xpost(1),xpost(2),'p')
hold off
title(sprintf('truestate=[%4.3f,%4.3f]',truestate(1),truestate(2)))

end

% fprintf('Error in init: %6.4f\n',norm(xinit - truestate))
% fprintf('Error in obs : %6.4f\n',norm(yobs  - truestate))
% fprintf('Error in EnKF: %6.4f\n',norm(xpost - truestate))
% fprintf('Error in PF  : %6.4f\n',norm(ppost - truestate))
% fprintf('Error in PFxR: %6.4f\n',norm(ppost1 - truestate))

errorarray(i_inflmu,i_n_ens,ix,:) = [...
    norm(xinit - truestate),...
    norm(yobs  - truestate),...
    norm(xpost - truestate),...
    norm(ppost - truestate),...
    norm(ppost1 - truestate)];

end
end
end

sound(y(1:18000),Fs)

%% A simple 3-dimensional case


for ix = 1:3
for iy = 1:3

ndim = 3;
dobs = 1;

Pf = 1.5*eye(ndim);
truestate = xlong(1:3,127*ix*iy);
xinit = mvnrnd(truestate,Pf)';

sigmaobs = 0.16;
Robs = sigmaobs.^2.* eye(dobs);

H = [1,1,1];
% H = eye(2);

yobs = mvnrnd(H * truestate, Robs,1)';
% disp(yobs)
% s = sum(yobs);
% yobs(1) = yobs(1) * (0.5 + 0.01*randn);
% yobs(2) = s - yobs(1);
% disp(yobs)

K = Pf * H' / (H * Pf * H' + Robs);

xpost = xinit + K * (yobs - H * xinit);
Pa = (eye(ndim) - K * H) * Pf;
%Pa = (Pa + Pa')/2;

nens = 1e7;

ens = mvnrnd(xpost,Pa,nens);

fprintf('H(x)=%4.4f, H(x_f)=%4.4f, yobs=%4.4f, H(x_a)=%4.4f\n',...
    H*truestate,H*xinit,yobs,H*xpost)

figure
for fi = 1:3
subplot(1,3,fi)
if fi==1
    ind = [2,3];
elseif fi==2
    ind = [1,3];
elseif fi==3
    ind = [1,2];
end
nbins=[100 100];
[N,C]=hist3(ens(:,ind),nbins);
lN = log(N);
lN(lN==-inf)=0;
contourf(C{1},C{2},lN',8)
% contourf(C{1},C{2},N',10)
hold on
scatter(xinit(ind(1)),xinit(ind(2)),'o')
scatter(truestate(1),truestate(2),'*')
scatter(xpost(ind(1)),xpost(ind(2)),'p')
hold off
end
end
end

%% A 2-pairwise 40-dimensional case


for ix = 1:2
    
for iy = 1:2

figure
    
ndim = 20;
dobs = 10;

Pf = 1.5*eye(ndim);
truestate = [1.6944    2.8240    2.0793   -3.2302    2.3776 ...
             8.0174    3.6318   -0.4771   -0.9462    2.2906 ...
             4.1460    8.7287    6.5204   -1.1531    6.2053 ...
             6.7588    1.3073    4.8381    1.4212   -3.8420]';
xinit = mvnrnd(truestate,Pf)';

sigmaobs = 0.16;
Robs = sigmaobs.^2.* eye(dobs);

% H = [1,1];
% H = zeros(dobs,ndim);
% for i=1:dobs
%     H(i,2*i-1:2*i) = 1;
% end

H = [eye(10),zeros(10,10)];

yobs = mvnrnd(H * truestate, Robs,1)';

for i=1:dobs/2
s = sum(yobs(2*i-1:2*i));
yobs(2*i-1) = yobs(2*i-1) * (0.5 + 0.1*randn);
yobs(2*i) = s - yobs(2*i-1);
end

K = Pf * H' / (H * Pf * H' + Robs);

xpost = xinit + K * (yobs - H * xinit);
Pa = (eye(ndim) - K * H) * Pf;
%Pa = (Pa + Pa')/2;

nens = 1e7;

ens = mvnrnd(xpost,Pa,nens);

for fi=1:5
subplot(2,5,fi)
nbins=[100 100];
[N,C]=hist3(ens(:,2*fi-1:2*fi),nbins);
lN = log(N);
lN(lN==-inf)=0;
contourf(C{1},C{2},lN',8)
hold on
scatter(xinit(2*fi-1),xinit(2*fi),'o')
scatter(truestate(2*fi-1),truestate(2*fi),'*')
scatter(xpost(2*fi-1),xpost(2*fi),'p')
hold off
end



end
end
%%
function y = two_pair_stupid(x)

[n_dim, n_col] = size(x);

y = zeros(n_dim,n_col);

for j=1:n_col
    for i=1:n_dim/2
        
        s = sum(x(2*i-1:2*i,j));
        y(2*i-1,j) = x(2*i-1,j) * (0.5+0.1*randn);
        y(2*i,  j) = s - y(2*i-1,j);
        
    end
end

end