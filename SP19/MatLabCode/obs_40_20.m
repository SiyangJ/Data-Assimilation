function y = obs_40_20(x)
%OBS Summary of this function goes here
%   x: each column is a state vector.
%   H: R^40 -> R^10 -> R^20
%   Projection followed by nonlinear transformation

n_col = size(x,2);

% projection step R^40 -> R^10

proj_ind = 1:10;
x_proj = x(proj_ind,:);

% nonlinear transformation R^10 -> R^20

y = zeros(20,n_col);

% TODO: Optimization in terms of efficiency.
for j=1:n_col
for i=1:5
        
    p = x_proj(2*i-1,j);
    q = x_proj(2*i,j);
    
    p2 = p^2;
    q2 = q^2;
    
    y(4*i-3:4*i,j) = [p*q;p2+q2;p2;q2];   
end
end

end

