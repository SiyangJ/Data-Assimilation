function dy = diff_eqs(t,y,Fc)
dy=zeros(2,1);
Fs=100;%*(sin(12*pi*t)).^2;
Fo=85; 
%Fc=100;%normrnd(100,100);
Ft=2.8; 
Fb=2;
cw=1;
Hml=6.3;
K=100;
aml=0.2;
L=9.5;
hc=10;
a=0.8;
b=0.6;
c=0.6;
d=0.2;
Ka=100;
Kb=100;
Kc=100;
Kd=100;
%dy(1) = (1-((aml+y(2))/2+tanh(y(1)./(L*hc)).*(aml-y(2))./2))*Fs-Fo+Fc-Ft*y(1)/(cw*Hml)+Fb;

dy(1)= (1-((.2+y(2))/2+((.2-y(2))/2).*tanh(y(1)/(9.5*.5)))).*((1205/12)+(1/2).*(165+50.*sqrt(3)).*sin(2.*pi.*t)/pi+(1/2).*(-237.*sqrt(3)-476).*cos(2.*pi.*t)/pi-(447/4).*sin(4.*pi.*t)/pi+(275/4).*sqrt(3).*cos(4.*pi.*t)/pi+35.*sin(6.*pi.*t)/pi+(73/3).*cos(6.*pi.*t)/pi-(39/8).*sin(8.*pi.*t)/pi-(55/8).*sqrt(3).*cos(8.*pi.*t)/pi)-((253/3)+(1/2).*(44.*sqrt(3)+64).*sin(2.*pi.*t)/pi+(1/2).*(116+40.*sqrt(3)).*cos(2.*pi.*t)/pi+20.*sin(4.*pi.*t)/pi+sqrt(3).*cos(4.*pi.*t)/pi-(52/3).*sin(6.*pi.*t)/pi-(16/3).*cos(6.*pi.*t)/pi-12.*sin(8.*pi.*t)/pi+(13/2).*sqrt(3).*cos(8.*pi.*t)/pi)+Fc-(337/120+(1/20).*(6.*sqrt(3)+5).*sin(2.*pi.*t)/pi+(1/20).*(12+3.*sqrt(3)).*cos(2.*pi.*t)/pi+(11/40).*sin(4.*pi.*t)/pi+(1/40).*sqrt(3).*cos(4.*pi.*t)/pi-(1/6).*sin(6.*pi.*t)/pi-(1/10).*cos(6.*pi.*t)/pi-(9/80).*sin(8.*pi.*t)/pi+(7/80).*sqrt(3).*cos(8.*pi.*t)/pi).*y(1)/(6.3)+2;

if( ( (aml+y(2))/2+tanh(y(1)/(L*hc)).*(aml-y(2))/2 )-0.6>0 ) 

dy(2) = (y(1).^2/(Ka^2)).*y(2).*(1-y(2)./a)+(Kb^2/(1+y(1).^2)).*y(2).*(1-y(2)./b); %+normrnd(0,0.5);

%elseif (((aml+y(2))/2+tanh(y(1)/(L*hc)).*(aml-y(2))/2 )-0.6==0)

%dy(2)=((y(1).^2/(Ka^2)).*y(2).*(1-y(2)./a)+(Kb^2/(1+y(1).^2)).*y(2).*(1-y(2)./b) ).*((0.5*(aml-y(2))*(sech(y(1)/(L*hc))).^2).* (1-((aml+y(2))/2+tanh(y(1)./(L*hc)).*(aml-y(2))./2))*Fs-Fo+Fc-Ft*y(1)/(cw*Hml)+Fb + (1-0.5*tanh(y(1)/(L*hc))).*( (Kc^2/(1+y(1).^2)).*y(2).*(1-y(2)./c)+(y(1).^2/Kd^2).*y(2).*(1-y(2)./d)))./( (1-0.5*tanh(y(1)/(L*hc))).*(  (Kc^2/(1+y(1).^2)).*y(2).*(1-y(2)./c)+(y(1).^2/Kd^2).*y(2).*(1-y(2)./d)-((y(1).^2/(Ka^2)).*y(2).*(1-y(2)./a)+(Kb^2/(1+y(1).^2)).*y(2).*(1-y(2)./b) ) )  ) + (1-(((0.5*(aml-y(2))*(sech(y(1)/(L*hc))).^2).* (1-((aml+y(2))/2+tanh(y(1)./(L*hc)).*(aml-y(2))./2))*Fs-Fo+Fc-Ft*y(1)/(cw*Hml)+Fb + (1-0.5*tanh(y(1)/(L*hc))).*( (Kc^2/(1+y(1).^2)).*y(2).*(1-y(2)./c)+(y(1).^2/Kd^2).*y(2).*(1-y(2)./d)))./( (1-0.5*tanh(y(1)/(L*hc))).*(  (Kc^2/(1+y(1).^2)).*y(2).*(1-y(2)./c)+(y(1).^2/Kd^2).*y(2).*(1-y(2)./d)-((y(1).^2/(Ka^2)).*y(2).*(1-y(2)./a)+(Kb^2/(1+y(1).^2)).*y(2).*(1-y(2)./b) ) ) ))).*((Kc^2/(1+y(1).^2)).*y(2).*(1-y(2)./c)+(y(1).^2/Kd^2).*y(2).*(1-y(2)./d));

else
dy(2) = (Kc^2/(1+y(1).^2)).*y(2).*(1-y(2)./c)+(y(1).^2/Kd^2).*y(2).*(1-y(2)./d); %+normrnd(0,0.5);
end
