function xout = lorenz96(~,x)
xout = ([x(2:end,:); x(1,:)] - [x(end-1:end,:); x(1:end-2,:)]).*[x(end,:); x(1:end-1,:)] - x - 16;
