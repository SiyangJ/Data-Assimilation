function x = stupid_inverse_20_10(y)
%STUPID_INVERSE_20_10 Summary of this function goes here
%   Detailed explanation goes here

n_col = size(y,2);

% nonlinear transformation R^10 -> R^20

x = zeros(10,n_col);

% TODO: Optimization in terms of efficiency.
for j=1:n_col
for i=1:5
        
    p = y(4*i-1,j);
    q = y(4*i,j);
    
    x(2*i-1:2*i,j) = [sqrt(p);sqrt(q)];
end
end


end

