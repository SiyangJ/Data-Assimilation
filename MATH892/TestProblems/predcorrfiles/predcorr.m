function [x,a,Del]=predcorr(x,a,Del,tol)
J = fxn(x,a);
tan = -J\fan(x,a);
a1=a+Del; x1=x+tan*Del; k=0; kmax=7;
f=fn(x1,a1); J=fxn(x1,a1); 
while (k<=kmax) && (norm(f)> tol) % stationary Newton
    s = -J\f;
    x1=x1+s; f=fn(x1,a1); k=k+1;
end
Del=2^((4-k)/3)*Del; % next step to be tried
if (norm(f) >tol)
    a1=a; 
    x1=x; % corrector failed, reset
else
x=x1; a=a1; 
end
