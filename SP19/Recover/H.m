function [hval,isterminal,direction]=H(t,y)
aml=0.2 ;
L=9.5;
hc=10;
hval=((aml+y(2))/2+tanh(y(1)/(L*hc)).*(aml-y(2))/2)-0.6;
isterminal=0;
direction=0;

