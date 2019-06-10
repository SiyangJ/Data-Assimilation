function [xdot] = simplexsq(~,x)
xdot = [x(2:end,:); x(1,:)];
    