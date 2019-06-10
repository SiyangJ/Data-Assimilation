function y = EW09M(t,x,Fc,deltaobs)
%EW09M Summary of this function goes here
%   Detailed explanation goes here
[~,Y,~,~,~,~] = disode45(@(s,y) diff_eqs(s,y,Fc) ,@H, [t,t+deltaobs],x);

y = Y(end,:)';

end

