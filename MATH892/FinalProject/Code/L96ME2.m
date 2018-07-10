function dx = L96ME2(t,x,F,xi,zeta,gamma)

dx = L96(t,x+xi,F)+zeta-gamma*(x.^2);

end