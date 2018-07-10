% Time independent LSWE with order 1 inertia for forecasting in EnKF

function dw = ens_lsw_ti_order1(~,w)

u_0 = 1;
epsilon = 0.01;
l = length(w);
dw = zeros(l,1);

for i=1:l/4
    dw(4*i-3) = -2*pi*u_0*sin(2*pi*w(4*i-3))*cos(2*pi*w(4*i-2))-epsilon*4*pi^2*u_0^2*sin(2*pi*w(4*i-3))*cos(2*pi*w(4*i-2));
    dw(4*i-2) = 2*pi*u_0*cos(2*pi*w(4*i-3))*sin(2*pi*w(4*i-2))-epsilon*4*pi^2*u_0^2*sin(2*pi*w(4*i-3))*cos(2*pi*w(4*i-2));
    dw(4*i-1) = 0;
    dw(4*i) = 0;
end