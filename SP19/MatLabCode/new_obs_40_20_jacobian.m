function H = obs_40_20_jacobian(x)

H = zeros(20,40);

for i=1:5
        
    p = x(2*i-1);
    q = x(2*i);
    
    %p2 = p^2;
    %q2 = q^2;
    
    H(4*i-3:4*i,2*i-1:2*i) = ...
    [2*(p+q),2*(p+q);...
     1      ,      1;...
     2*(p+q),2*(p+q);...
     1      ,      1];
    
    %y(4*i-3:4*i) = [p*q;p2+q2;p2;q2];   
end

end