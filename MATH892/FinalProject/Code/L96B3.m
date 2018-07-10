function dx = L96B3(t,x,F,n)

dx = zeros(3*n,1);

% Evolve the bias estimators
% For simplicity, just constant
dx(n+1:end) = x(n+1:end);

% Evolve the system estimator
dx(1:n) = L96(t,x(1:n),F)+dx(n+1:2*n);

end