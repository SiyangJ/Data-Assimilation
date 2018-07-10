n=20;

xs = zeros(1,n);
as = zeros(1,n);

xs(1) = 1/sqrt(2);
as(1) = 1/sqrt(2);

Del = 0.1;

tol = 1e-5;

for k=2:n
    [xs(k),as(k),Del]=predcorr(xs(k-1),as(k-1),Del,tol);
end

%%
figure
hold on
plot(xs,'-o')
plot(as,'-o')
hold off
legend('x','a')

%%
n=100;

xs = zeros(1,n);
as = zeros(1,n);

xs(1) = 1/sqrt(2);
as(1) = 1/sqrt(2);

Del = 0.1;

tol = 1e-8;

for k=2:n
    [xs(k),as(k),Del]=pseudoarc(xs(k-1),as(k-1),Del,tol);
end

%%
figure
hold on
plot(xs,'-o')
plot(as,'-o')
hold off
legend('x','a')