hc=10;
a=@(E,AM) (0.2+AM)/2+tanh(E/(9.5*hc)).* (0.2-AM)/2 - 0.6;

[X,Y] = meshgrid(-200:0.1:250,0.2:1e-4:0.8);
A = a(X,Y);

tol = 1e-4;
A0 = abs(A)<tol;
X0 = X(A0);
Y0 = Y(A0);
%plot(X0,Y0)
fprintf('Finish part 1')
%%

%options = disodeset('RelTol',1e-2,'AbsTol',1e-3);

Fc=60;
x0=[-200; 0.7];
%[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, [0 50],x0);
[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, [0 50],x0);%, options);
hc=10;
a=((0.2+y(:,2))/2+tanh(y(:,1)/(9.5*hc)).*(0.2-y(:,2))/2);
Ci=1-((0.8-(.5*y(:,2)+.5*a))./0.6);%((y(:,2)-0.2)./(y(:,2))).*(a./y(:,2))+.2*(y(:,2)-0.2)./y(:,2);%(a-0.2)./y(:,2)+.2*(y(:,2)-0.2)./y(:,2);%(0.5+0.5.*tanh(-y(:,1)./10)).*a./y(:,2);%
Cp=(1-a./y(:,2));%Ci.*max(0,(0.6-a)./0.6);%(0.5+0.5*tanh(0.6-a./0.6))%(1-a./(y(:,2)));
Csat= max(0,Ci-Cp);
[m,~]=size(a);
Rad(:,1)=abs(y(:,1).*y(:,2));
Rad(:,2)=y(:,2)-a;
Rad(:,3)=(a).*abs(y(:,1));
Rad(:,4)=(.5+.4.*tanh((-(y(:,1)-50)/10))).*((y(:,1))+273.15);
Rad(:,5)=Cp.*Ci;
figure('NumberTitle', 'off', 'Name', sprintf('Fc = %d',i),...
       'rend','painters','pos',[10 10 900 600])
subplot(2,3,1)
plot(y(:,1),y(:,2))
plot(X0,Y0)
title('State space trajectory')
subplot(2,3,2)
plot(t,y(:,1))
title('Energy versus time')
subplot(2,3,3)
plot(y(:,1),y(:,2))
hold on
plot(y(:,1),a);
title('$$\alpha_m$$ and $$\alpha$$ versus $$E$$','interpreter','latex')
legend({'\alpha_m$$','$$\alpha'})
subplot(2,3,4)
plot(y(:,1),Ci)
title('Ice concentration versus energy')
subplot(2,3,5)
plot(y(:,1),Cp)
title('Pond fraction versus energy')
subplot(2,3,6)
plot(y(:,1),Ci-Csat)
title('Ice concentration error versus energy')
