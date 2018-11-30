import numpy as np
import math
import config

def ObsVer1(X):
    hc = 10 #TODO
    a =((0.2+X(:,2))/2+tanh(X[:,0]/(9.5*hc)).*(0.2-X(:,2))/2);
    Ci=1-((0.8-(.5*X(:,2)+.5*a))./0.6);
    Cp=(1-a./X(:,2));
    Csat= max(0,Ci-Cp);
    [m,~]=size(a);
    Rad = np.zeros(,5)
    Rad[:,0]=abs(X[:,0].*X(:,2));
    Rad(:,2)=X(:,2)-a;
    Rad(:,3)=(a).*abs(X[:,0]);
    Rad(:,4)=(.5+.4.*tanh((-(X[:,0]-50)/10))).*((X[:,0])+273.15);
    Rad(:,5)=Cp.*Ci;
    return Rad