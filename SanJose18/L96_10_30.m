function dx = L96_10_30(~,x,F)
% The first 10 variables and the 30 latter variables evolve independently.
% Experiment to see if uncorrelated variables make
% h_10 better than h_40
% The equations are similar to Lorenz 96

dx = zeros(length(x),1);

n = 10;

% dx(n) = (x(1)-x(n-2))*x(n-1)-x(n)+F;
% dx(1) = (x(2)-x(n-1))*x(n)-x(1)+F;
% dx(2) = (x(3)-x(n))*x(1)-x(2)+F;
% 
% dx(3:n-1) = (x(4:n)-x(1:n-3)).*x(2:n-2)-x(3:n-1)+F;
% 
% dx(end) = (x(n+1)-x(end-2))*x(end-1)-x(end)+F;
% dx(n+1) = (x(n+2)-x(end-1))*x(end)-x(n+1)+F;
% dx(n+2) = (x(n+3)-x(end))*x(n+1)-x(n+2)+F;
% 
% dx(n+3:end-1) = (x(n+4:end)-x(n+1:end-3)).*x(n+2:end-2)-x(n+3:end-1)+F;

dx(1:n) = L96(0,x(1:n),F);
dx(n+1:end) = L96(0,x(n+1:end),100);


end

