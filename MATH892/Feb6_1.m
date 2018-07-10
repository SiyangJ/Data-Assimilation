Dim = [5,10,20,50,80,100,200];
d = length(Dim);
errs = zeros(d,4);
time = zeros(d,5);

tol = 1e-7;

for i = 1:d
n = Dim(i);
Q = orth(randn(n));
D = diag(abs(randn(n,1))+0.5);
A = Q*D*Q';

x = randn(n,1);

b = A*x;
tic
R = chol(A);
time(i,5) = toc;
tic
x_ = A\b;
time(i,1) = toc;
tic
Rx = R'\b;
x_c = R\Rx;
time(i,2) = toc;

Af = @(x)A*x;
tic
x_gf = pcg(Af,b,tol);
time(i,4) = toc;
tic
x_g = pcg(A,b,tol);
time(i,3) = toc;

err_ = norm(x_-x);
err_c = norm(x_c-x);
err_g = norm(x-x_g);
err_gf = norm(x-x_gf);

errs(i,:) = [err_,err_c,err_g,err_gf];

end
%%
loglog(Dim,time)
legend('slash','chol','pcg','pcg-f','chol-fac')
%%
loglog(Dim,errs)
legend('slash','chol','pcg','pcg-f')