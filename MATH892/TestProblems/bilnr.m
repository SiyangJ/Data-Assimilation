function b = bilnr(n,u,v,par)
a = par(1);
uu(1:n) = u(1:n);
uu(n+1:3*n) = 0.0;

for k=1:2*n
    sum = 0.0;
    for j=1:n
        if k > j
            kjsgn=1;
        elseif k < j
            kjsgn=-1;
        else
            kjsgn=0;
        end
        if k==j
            sum = sum + j*v(j)*(uu(j+k));
        else
            sum = sum + j*v(j)*(uu(j+k)+kjsgn*uu(abs(k-j)));
        end
    end
    b(k) = a*sum/2.0;
end
