figure('NumberTitle', 'off', 'Name', 'A bunch of trajectories',...
       'rend','painters','pos',[10 10 900 600])
hold on   
irange = 1:500;
labels = cell(1,length(irange));
tspan = [0,100];
for i=irange

Fc=0;
x0=[min(-200 + randn()*75,0.01); min(max(0.75 + randn()*0.4,0.01),0.99)];
[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, tspan,x0);
% hc=10;
% a=((0.2+y(:,2))/2+tanh(y(:,1)/(9.5*hc)).*(0.2-y(:,2))/2);
% Ci=1-((0.8-(.5*y(:,2)+.5*a))./0.6);
% Cp=(1-a./y(:,2));
% Csat= max(0,Ci-Cp);
% [m,~]=size(a);
% Rad(:,1)=abs(y(:,1).*y(:,2));
% Rad(:,2)=y(:,2)-a;
% Rad(:,3)=(a).*abs(y(:,1));
% Rad(:,4)=(.5+.4.*tanh((-(y(:,1)-50)/10))).*((y(:,1))+273.15);
% Rad(:,5)=Cp.*Ci;

plot(y(:,1),y(:,2))
labels(i) = cellstr(sprintf('x0=[%4.2f;%4.4f]',x0(1),x0(2)));

end
hold off
%legend(labels)
%%
hc=10;
a=@(E,AM) (0.2+AM)/2+tanh(E/(9.5*hc)).* (0.2-AM)/2 - 0.6;

[X,Y] = meshgrid(-200:1:250,0.2:1e-3:0.8);
A = a(X,Y);

tol = 2e-4;
A0 = abs(A)<tol;
X0 = X(A0);
Y0 = Y(A0);
plot(X0,Y0)
%%

[XX,YY] = meshgrid(-200:15:250,0.2:2e-2:0.8);

figure('NumberTitle', 'off', 'Name', 'A bunch of trajectories',...
       'rend','painters','pos',[10 10 900 600])
hold on   
irange = 1:size(XX,1);
jrange = 1:size(XX,2);
tspan = [0,100];
for i=irange
for j=jrange

Fc=0;
x0=[XX(i,j);YY(i)];
[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, tspan,x0);
% hc=10;
% a=((0.2+y(:,2))/2+tanh(y(:,1)/(9.5*hc)).*(0.2-y(:,2))/2);
% Ci=1-((0.8-(.5*y(:,2)+.5*a))./0.6);
% Cp=(1-a./y(:,2));
% Csat= max(0,Ci-Cp);
% [m,~]=size(a);
% Rad(:,1)=abs(y(:,1).*y(:,2));
% Rad(:,2)=y(:,2)-a;
% Rad(:,3)=(a).*abs(y(:,1));
% Rad(:,4)=(.5+.4.*tanh((-(y(:,1)-50)/10))).*((y(:,1))+273.15);
% Rad(:,5)=Cp.*Ci;

plot(y(:,1),y(:,2))

end
end
plot(X0,Y0)
hold off

%%
x0=[-243.87249;0.7479624];
tspan = [0,50];
Fc = 0;
[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, tspan,x0);

figure
plot(y(:,1),y(:,2))