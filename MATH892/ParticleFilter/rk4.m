function [x] = rk4(model,Tval,dt,IC)
%attempt to integrate between the two entries of Tval with a time step of
%dt. 
nsteps = ceil((Tval(2)-Tval(1))/dt);
hstep = (Tval(2)-Tval(1))/nsteps; %numerical time step
x = IC;
% x = zeros(nsteps,length(IC)); %initialise solution vector
% x(1,:) = IC;
% xtemp = IC; %xtemp holds only the state at the current time - to speed computation

for n=1 : nsteps-1
    t = Tval(1)+(n-1)*hstep; %current time    
    k1 = hstep*feval(model,t,x); %rk4 uses 4 `increments' labelled k1 to k4
    k2 = hstep*feval(model,t,x+k1/2);
    k3 = hstep*feval(model,t,x+k2/2);
    k4 = hstep*feval(model,t,x+k3);
    
    x = x + (k1+2*k2+2*k3+k4)/6; %rk4 update for state
    %x(n+1,:) = xtemp; %store current state
end