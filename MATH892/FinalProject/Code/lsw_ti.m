% Time-independent LSWE without inertia

function dz = lsw_ti(~,z)

u_0 = 1;
dz = zeros(4,1);

dz(1) = -2*pi*u_0*sin(2*pi*z(1)).*cos(2*pi*z(2));       % z(1) = x
dz(2) = 2*pi*u_0*cos(2*pi*z(1)).*sin(2*pi*z(2));        % z(2) = y
dz(3) = 0;                                              % z(3) = U
dz(4) = 0;                                              % z(4) = V