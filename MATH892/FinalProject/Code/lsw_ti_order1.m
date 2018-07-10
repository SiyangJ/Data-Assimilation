% Time independent LSWE with order 1 inertia for forecasting in EnKF

function dw = lsw_ti_order1(~,w)

u_0 = 1;
epsilon = 0.01;
dw = zeros(4,1);

    dw(1) = -2*pi*u_0*sin(2*pi*w(1))*cos(2*pi*w(2))-epsilon*4*pi^2*u_0^2*sin(2*pi*w(1))*cos(2*pi*w(2));
    dw(2) = 2*pi*u_0*cos(2*pi*w(1))*sin(2*pi*w(2))-epsilon*4*pi^2*u_0^2*sin(2*pi*w(1))*cos(2*pi*w(2));
    dw(3) = 0;
    dw(4) = 0;
end