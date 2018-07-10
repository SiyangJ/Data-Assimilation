function fu = kse(t, u, par)
%RHS of KSE
a=par(1); %a=32 typical value.
n = length(u);
uu(1:n) = u(1:n);
uu(1+n:3*n) = zeros;
b = bilnr(n,uu,uu,par);
for j=1:n
    fu(j) = -4*j^4*u(j)+a*j^2*u(j) - b(j);
end

end
