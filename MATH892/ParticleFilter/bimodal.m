function xdot = bimodal(~,x)
[p1,p2]=size(x);
xdot = -x.*(x-1).*(x+1);