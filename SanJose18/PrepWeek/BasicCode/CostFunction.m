%% Objective function definition
function out=CostFunction(u,y,gamma,alpha,m0,C0,J)

Phi=1/2/C0*(u-m0)^2;
uu(1)=u;
for j=1:J
    uu(j+1) = alpha*sin(uu(j));
    Phi=Phi+1/2/gamma^2*(y(j)-uu(j+1))^2;
end
out=Phi;
end



