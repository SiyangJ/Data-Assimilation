function [x1,P1] = Predict(A,Q,x,P)
%PREDICT Summary of this function goes here
%   Detailed explanation goes here

x1 = A*x;
P1 = A*P*A'+Q;

end

