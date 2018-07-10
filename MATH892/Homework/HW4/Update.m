function [x1,P1] = Update(H,R,x,P,y)
%UPDATE Summary of this function goes here
%   Detailed explanation goes here

K = P*H'*(H*P*H'+R)^-1;
x1 = x+K*(y-H*x);
P1 = P-K*H*P;

end

