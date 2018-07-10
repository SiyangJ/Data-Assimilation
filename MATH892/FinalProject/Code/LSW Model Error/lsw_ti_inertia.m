% Time independent LSWE with inertia

function dw = lsw_ti_inertia(~,w)

u_0 = 1;
epsilon = 0.01;
dw = zeros(4,1);

dw(1) = epsilon*w(3);                                       % w(1) = x
dw(2) = epsilon*w(4);                                       % w(2) = y
dw(3) = -2*pi*u_0*sin(2*pi*w(1))*cos(2*pi*w(2))-w(3);       % w(3) = U
dw(4) = +2*pi*u_0*cos(2*pi*w(1))*sin(2*pi*w(2))-w(4);       % w(4) = V