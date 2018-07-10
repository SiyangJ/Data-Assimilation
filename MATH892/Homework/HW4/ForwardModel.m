function next = ForwardModel(A,Q,cur)
%FORWARDMODEL Summary of this function goes here
%   Detailed explanation goes here

next = A*cur+mvnrnd(zeros(1,length(cur)),Q)';

end

