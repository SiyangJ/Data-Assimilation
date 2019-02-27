
x0=[-250; 0.75];

Fc=0;
for i=0
Fc=i;
x0=[-200; 0.75];
%[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, [0 50],x0);
[t,y,tdis,ydis,idis,stats]=disode45(@(t,y) diff_eqs(t,y,Fc) ,@H, [0 50],x0);
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
title('State space trajectory')
subplot(2,3,2)
plot(t,y(:,1))
title('Energy versus time')
subplot(2,3,3)
plot(y(:,1),y(:,2))
hold on
plot(y(:,1),a);
title('$$\alpha_m$$ and $$\alpha$$ versus $$E$$','interpreter','latex')
subplot(2,3,4)
plot(y(:,1),Ci)
title('Ice concentration versus energy')
subplot(2,3,5)
plot(y(:,1),Cp)
title('Pond fraction versus energy')
subplot(2,3,6)
plot(y(:,1),Ci-Csat)
title('Ice concentration error versus energy')
%figure;
%plot(ydis(:,1),ydis(:,2))

%pause
%x0=[y(i,1),y(i,2)]
state(:,1)=y(:,1);
state(:,2)=y(:,2);
state(:,3)=a;
state(:,4)=Ci;
state(:,5)=Cp;

% Rade(:,1)=Rad(:,1)+normrnd(0,mean(Rad(:,1))/max(Rad(:,1)),m,1);
% Rade(:,2)=Rad(:,2)+normrnd(0,mean(Rad(:,2))/max(Rad(:,2)),m,1);
% Rade(:,3)=Rad(:,3)+normrnd(0,mean(Rad(:,3))/max(Rad(:,3)),m,1);
% Rade(:,4)=Rad(:,4)+normrnd(0,mean(Rad(:,4))/max(Rad(:,4)),m,1);
% Rade(:,5)=Rad(:,5)+normrnd(0,mean(Rad(:,5))/max(Rad(:,5)),m,1);

% Rade(:,1)=Rad(:,1)+normrnd(0,mean(Rad(:,1))/2,m,1);
% Rade(:,2)=Rad(:,2)+normrnd(0,mean(Rad(:,2))/2,m,1);
% Rade(:,3)=Rad(:,3)+normrnd(0,mean(Rad(:,3))/2,m,1);
% Rade(:,4)=Rad(:,4)+normrnd(0,mean(Rad(:,4))/2,m,1);
% Rade(:,5)=Rad(:,5)+normrnd(0,mean(Rad(:,5))/2,m,1);
% 
% 
% dlmwrite(strcat('Data/statevariables',num2str(Fc),'.txt'),state)
% dlmwrite(strcat('Data/radiances',num2str(Fc),'.txt'),Rad)
% dlmwrite(strcat('Data/radiances_error',num2str(Fc),'.txt'),Rade)

clear
end;

