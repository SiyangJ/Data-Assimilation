% Time independent LSWE without inertia for forecasting in EnKF

function dw = ens_lsw_ti(~,w)

u_0 = 1;
l = length(w);
dw = zeros(l,1);

for i=1:l/4
    dw(4*i-3) = -2*pi*u_0*sin(2*pi*w(4*i-3))*cos(2*pi*w(4*i-2));       % w(4*i-3) = x
    dw(4*i-2) = 2*pi*u_0*cos(2*pi*w(4*i-3))*sin(2*pi*w(4*i-2));        % w(4*i-2) = y
    dw(4*i-1) = 0;                                              % w(4*i-1) = U
    dw(4*i) = 0;                                                % w(4*i) = V
end