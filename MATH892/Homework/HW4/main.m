%% Setting up constants
dt = 0.3;
A = eye(4);
A(1,1)=1;
A(2,2)=1;
A(3,3)=1;
A(4,4)=1;
A(1,3)=dt;
A(2,4)=dt;
H = [1,0,0,0;0,1,0,0];
q1c = 0.03;
q2c = 0.03;
Q = [q1c*dt^3/3,0,q1c*dt^2/2,0;...
    0,q2c*dt^3/3,0,q2c*dt^2/2;...
    q1c*dt^2/2,0,q1c*dt,0;...
    0,q2c*dt^2/2,0,q2c*dt];
s1 = 0.02;
s2 = 0.02;
R = [s1,0;0,s2];


x0a = [0;0;1;1];
P0a = eye(4);

%% Generating trues
rng(1);
x0t = [0.5;0.3;-1;-0.2];
N1 = 50;
N2 = 40;
xt = zeros(4,N1+1);
xt(:,1)=x0t;
for i=2:N1+1
    xt(:,i) = ForwardModel(A,Q,xt(:,i-1));
end
yt = zeros(2,N2+1);
for i=1:N2+1
    yt(:,i) = Observe(H,R,xt(:,i));
end

%% Pure Prediction
rng(1);
xp = zeros(4,N1+1);
Pp = zeros(4,4,N1+1);
xp(:,1)=x0a;
Pp(:,:,1)=P0a;
for i=2:N1+1
    [xp(:,i),Pp(:,:,i)]=Predict(A,Q,xp(:,i-1),Pp(:,:,i-1));
end

errp = abs(xp-xt);

%% Analysis
rng(1);
xf = zeros(4,N1+1);
Pf = zeros(4,4,N1+1);
xa = zeros(4,N1+1);
Pa = zeros(4,4,N1+1);
xa(:,1)=x0a;
Pa(:,:,1)=P0a;
for i=2:N2+1
    [xf(:,i),Pf(:,:,i)] = Predict(A,Q,xa(:,i-1),Pa(:,:,i-1));
    [xa(:,i),Pa(:,:,i)] = Update(H,R,xf(:,i),Pf(:,:,i),yt(:,i));
end

%% Extrapolating
i=N2+2;
[xf(:,i),Pf(:,:,i)] = Predict(A,Q,xa(:,i-1),Pa(:,:,i-1));
for i=N2+3:N1+1
    [xf(:,i),Pf(:,:,i)] = Predict(A,Q,xf(:,i-1),Pf(:,:,i-1));
end
%% Plot
for i=1:4
    subplot(2,2,i)
    hold on
    plot(0:N1,xt(i,:))
    errorbar(0:N1,xf(i,:),squeeze(Pf(i,i,:)))
    hold off
    legend('true','analysis')
    title(sprintf('Component %d',i))
end

