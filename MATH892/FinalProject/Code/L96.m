function dx = L96(~,x,F)

dx = zeros(length(x),1);

dx(end) = (x(1)-x(end-2))*x(end-1)-x(end)+F;
dx(1) = (x(2)-x(end-1))*x(end)-x(1)+F;
dx(2) = (x(3)-x(end))*x(1)-x(2)+F;

dx(3:end-1) = (x(4:end)-x(1:end-3)).*x(2:end-2)-x(3:end-1)+F;

end