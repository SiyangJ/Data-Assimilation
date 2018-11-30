function yout = Mdefault(x0,t0,t1,Fc)
%M Summary of this function goes here
%   Detailed explanation goes here
Fs=100;
Fo=85; 
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
Ka=150;
Kb=100;
Kc=50;
Kd=100;

[~,y,~,~,~,~]=disode45(@(t,y) diff_eqs_param(t,y,...
    Fc,Fs,Fo,Ft,Fb,...
    cw,Hml,K,aml,L,hc,...
    a,b,c,d,...
    Ka,Kb,Kc,Kd),...
    @H, [t0 t1],x0);

yout = y(end,:);
end
