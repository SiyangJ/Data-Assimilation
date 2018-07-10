% Plot inertial LSW trajectories

close all;

close all;

T_i = 0; T_f = 50;
h = 0.001;
comp_time = T_i:h:T_f;                        % computational time domain
comp_len = length(comp_time);
u_0 = 1;

ic1 = [0.; 0.2; 1; 1];
ic2 = [0.95; 0.05; 1; 1];
ic3 = [0.1; 0.9; 1; 1];
ic4 = [0.8; 0.7; 1; 1];

truth = zeros(16,comp_len);
truth(:,1) = [ic1; ic2; ic3; ic4];

figure
hold on

[x,y] = meshgrid(0:0.05:1,0:0.05:1);
u = -2*pi*u_0*sin(2*pi*x).*cos(2*pi*y);
v = 2*pi*u_0*cos(2*pi*x).*sin(2*pi*y);

for j=1:comp_len-1
    truth(:,j+1) = rk4('ens_lsw_ti',comp_time(j),comp_time(j+1),h,truth(:,j));
end
plot(truth(1,:),truth(2,:),'r')
%plot(truth(13,:),truth(14,:),'g')
title('Inertial Trajectories in the Cellular Flow Field')
xlabel('x')
ylabel('y')
quiver(x,y,u,v)
%leg = legend('ic = (0.1,0.2)','ic = (0.8, 0.7)');
%set(leg,'FontSize',18)