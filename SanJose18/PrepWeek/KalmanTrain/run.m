% Finish this .m file to get the K filter running!
load traindata.mat
% is=[0;100]; % the inital condition, we will take 0 to be where the train enters the forest going 100
% Q= ; % What should our process covariance matrix be?
% B0= ; % What guess might we take for our first error covariance matrix
% R= ;% What sould our observation coverianve matrix be
% H= ;% What would our observation operator be since we are directely observing both x and v?
% A= ;% What matrix will move our model forward.

is=[0;100];
Q=[(15/60)^2,225/60;225/60,225];
R=[100,0;0,100];
H=[1,0;0,1];
A=[1,1/60;0,1];
B0=eye(2,2);

[xk1,Bk1]=predict(is,B0,A,Q); % here we make our first prediction
y(1,1)=xGPS(1); % load the first observations
y(2,1)=vGPS(1); % load the first observations
[xa,Ba]=Kupdate(xk1,y,Bk1,H,R); % first analysis update
xk(1)=xk1(1); % store the first predictions
vk(1)=xk1(2); % store the first predictions
Xa(1)=xa(1); % store the first analysis
Va(1)=xa(2); % store the first analysis

for i=2:60 % this loop continues the processcl
    [xk1,Bk1]=predict(xa,Ba,A,Q);
    xk(i)=xk1(1);
    vk(i)=xk1(2);
    y(1,1)=xGPS(i);
    y(2,1)=vGPS(i);
    [xa,Ba]=Kupdate(xk1,y,Bk1,H,R);
    Xa(i)=xa(1);
    Va(i)=xa(2);
end

% Note for plot stuff you will have to hit enter after the inital plot to
% continue plotting, messages describing the plot will apper in the command
% window.
plotstuff;

 