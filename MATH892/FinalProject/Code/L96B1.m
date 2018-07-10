function dx = L96B1(t,x,F,n)

dx = zeros(2*n,1);

% Evolve the bias estimator
% For simplicity, just constant
dx(n+1:end) = 0;

% Evolve the system estimator
dx(1:n) = L96(t,x(1:n),F)+dx(n+1:end);

end