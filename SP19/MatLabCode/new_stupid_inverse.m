function y = new_stupid_inverse(x)
%NEW_OBS_40_20 Summary of this function goes here
%   x: each column is a state vector.
%   H: R^40 -> R^20
%   Projection followed by nonlinear transformation

n_col = size(x,2);

% projection step R^40 -> R^10

proj_ind = 1:10;
x_proj = x(proj_ind,:);

% nonlinear transformation R^10 -> R^20

y = zeros(10,n_col);

% TODO: Optimization in terms of efficiency.
for j=1:n_col
for i=1:5
        
    p = x_proj(2*i-1,j);
    q = x_proj(2*i,j);
    
    s = p+q;
    
    p = p * (0.5 + 0.1 * randn);
    %p = p + 3 * randn;
    q = s-p;
    
    % y(i,j) = sin(p^2 + q^2) + cos(p+q);
    y(2*i-1:2*i,j) = [p;q];
    
    % y(4*i-3:4*i,j) = [sin(p2+q2);cos(p2+q2);sin(p2);cos(p2)];
end
end

end

