function dx = L96(~,x,F)

% Summary:
% x' = L(x)
% Input: 
%        ~: Supposed to be time, since Lorenz 96' is not time dependent, we use ~ instead.
%        x: The state variable, x, a 40-dimensional vector.
%        F: The major parameter in Lorenz 96' model, also a 40-dimensional vector.
%           Usually, perple use F=[8;8;...;8], and it's not very chaotic.
%           I used sinusoidal wave to produce more complicated behavior,
%           suggested by Baek et al 2006, in the context of model error.
% Output:
%        dx: L(x), the derivative, a 40-dimensional vector

% Note:
%        Numerically, people usually use RK4 to forward the model.
%        We can code our own RK4 solver.
%        It's more effiencient, but less accurate.
%        But for the time being, 
%        I think we should probably just use the built-in solver ode45 instead.

dx = zeros(length(x),1);

dx(end) = (x(1)-x(end-2))*x(end-1)-x(end)+F;
dx(1) = (x(2)-x(end-1))*x(end)-x(1)+F;
dx(2) = (x(3)-x(end))*x(1)-x(2)+F;

dx(3:end-1) = (x(4:end)-x(1:end-3)).*x(2:end-2)-x(3:end-1)+F;

end