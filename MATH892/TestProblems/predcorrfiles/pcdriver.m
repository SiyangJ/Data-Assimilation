tol=1.E-8;
x=1/sqrt(2);
a=1/sqrt(2);
Del=0.1;
for k=1:20
[x,a,Del]=predcorr(x,a,Del,tol)
end
