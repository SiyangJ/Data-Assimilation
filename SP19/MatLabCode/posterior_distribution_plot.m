%% A simple 2-dimensional case

for ix = 1:1
    
for iy = 1:1
rng(19970215)
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

subplot(5,5,(ix-1)*5+iy)
nbins=[100 100];
[N,C]=hist3(ens,nbins);
%[N,C] =hist3w(ens,ones(nens,1)/nens,nbins);
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

%% A 2-dimensional case with ensembles


for inflmu = [0.05,0.1,0.2,0.4,0.8]
    
for n_ens = [10,20,100,1e4]
figure('NumberTitle', 'off', 'Name', sprintf('inflmu=%3.2f, nens=%d',inflmu,n_ens));
for ix = 1:6
ndim = 2;
dobs = 2;

Pf = 1.5*eye(ndim);
truestate = xlong(ix+10:ix+11,127*ix);
xinit = mvnrnd(truestate,Pf)';
xens = mvnrnd(xinit,Pf,n_ens)';

sigmaobs = 0.16;
Robs = sigmaobs.^2.* eye(dobs);

% H = [1,1];
H = eye(2);

h = @two_pair_stupid;

hidden_obs = H * truestate;
yobs = mvnrnd(h(truestate), Robs,1)';

pf = cov(xens');
pf = pf+inflmu*trace(pf)/ndim*eye(ndim);
pfht = pf*H';   
K = pfht*pinv(H*pfht+Robs);  

xenspost = xens + K * (mvnrnd(yobs,Robs,n_ens)' - H * xens);
xpost = mean(xenspost,2);
Pa = cov(xenspost');

nens = 1e7;

subplot(3,2,ix)
nbin = min(floor(n_ens/5),100);
nbins=[nbin,nbin];
[N,C]=hist3(xenspost',nbins);
lN = log(N);
lN(lN==-inf)=0;
nlevel = min(nbin,8);
contourf(C{1},C{2},lN',nlevel)
% contourf(C{1},C{2},N',10)
hold on
scatter(xinit(1),xinit(2),'o')
scatter(truestate(1),truestate(2),'*')
scatter(xpost(1),xpost(2),'p')
hold off
title(sprintf('truestate=[%4.3f,%4.3f]',truestate(1),truestate(2)))
end
end
end

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