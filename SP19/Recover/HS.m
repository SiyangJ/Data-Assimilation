function HS=HS(y,Fc)
aml=0.2;
L=9.5;
hc=10;
HS=((aml+y(:,2))/2+tanh(y(:,1)/(L*hc)).*(aml-y(:,2))/2)-0.6;