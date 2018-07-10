function [x,a,Del]=pseudoarc(x,a,Del,tol)
J = df(x,a);
v = null(J);

if (a>0 && v(1)>0) || (a<=0 && v(1)<0)
    v=-v;
end
a1=a+Del*v(2); 
x1=x+Del*v(1); 
k=0; kmax=7;
f=fn(x1,a1); 
while (k<=kmax) && (norm(f)> tol) % stationary Newton
    J(1,1)=fxn(x1,a1); 
    J(1,2)=fan(x1,a1);
    s = -pinv(J)*f;
    x1=x1+s(1); 
    a1=a1+s(2);
    f=fn(x1,a1); 
    k=k+1;
end
Del=2^((4-k)/3)*Del; % next step to be tried
if (norm(f) >tol)
    a1=a; 
    x1=x; % corrector failed, reset
else
    x=x1; 
    a=a1; 
end
