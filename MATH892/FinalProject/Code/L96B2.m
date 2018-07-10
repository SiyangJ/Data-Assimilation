function dx = L96B2(t,x,F,n)

dx = zeros(2*n,1);

% Evolve the bias estimators
% For simplicity, just constant
dx(n+1:end) = 0;

% Evolve the system estimator
dx(1:n) = L96(t,x(1:n),F);

end