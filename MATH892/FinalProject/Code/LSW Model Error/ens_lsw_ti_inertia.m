% Time independent LSWE with inertia for forecasting in EnKF

function dw = ens_lsw_ti_inertia(~,w)

u_0 = 1;
epsilon = 0.01;
l = length(w);
dw = zeros(l,1);

for i=1:l/4
    dw(4*i-3) = epsilon*w(4*i-1);
    dw(4*i-2) = epsilon*w(4*i);
    dw(4*i-1) = -w(4*i-1)-2*pi*u_0*sin(2*pi*w(4*i-3))*cos(2*pi*w(4*i-2));
    dw(4*i) = -w(4*i)+2*pi*u_0*cos(2*pi*w(4*i-3))*sin(2*pi*w(4*i-2));
end