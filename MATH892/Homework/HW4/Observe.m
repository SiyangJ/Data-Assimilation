function y = Observe(H,R,x)
%OBSERVE Summary of this function goes here
%   Detailed explanation goes here

y = H*x+mvnrnd(zeros(1,size(H,1)),R)';

end

