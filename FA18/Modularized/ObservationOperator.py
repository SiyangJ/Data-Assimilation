import numpy as np
import math
import config

def TrueObsVer1(X):
    hc = 10 #TODO
    a =((0.2+X[:,1])/2+np.tanh(X[:,0]/(9.5*hc)) * (0.2-X[:,1])/2)
    Ci = 1-((0.8-(0.5*X[:,1]+.5*a))/0.6)
    Cp = (1 - a / X[:,1])
    Csat = np.maximum(0,Ci-Cp)
    m = a.shape[0]
    Rad = np.zeros([m,5])
    Rad[:,0] = np.absolute(X[:,0] * X[:,1])
    Rad[:,1] = X[:,1] - a
    Rad[:,2] = a * np.absolute(X[:,0])
    Rad[:,3] = (0.5+0.4 *np.tanh((-(X[:,0]-50)/10))) * ((X[:,0])+273.15)
    Rad[:,4] = Cp * Ci
    return Rad

def StupidObsVer1(X):
    hc = 10 #TODO
    a =((0.2+X[:,1])/2+np.tanh(X[:,0]/(9.5*hc)) * (0.2-X[:,1])/2)
    Ci = 1-((0.8-(0.5*X[:,1]+.5*a))/0.6)
    Cp = (1 - a / X[:,1])
    Csat = np.maximum(0,Ci-Cp)
    return Csat

def EW09default(x0,t0,t1,Fc):
    import matlab as ml

    eng = config.MATLAB_ENGINE
    return eng.Mdefault(eng.transpose(ml.double(x0)),float(t0),float(t1),float(Fc))

def EW09(x0,t0,t1,Fc,Fs,Fo,Ft,Fb,
         cw,Hml,K,aml,L,hc,
         a,b,c,d,
         Ka,Kb,Kc,Kd):
    import matlab as ml

    eng = config.MATLAB_ENGINE
    return eng.Mdefault(eng.transpose(ml.double(x0)),float(t0),float(t1),float(Fc),
                       Fs,Fo,Ft,Fb,
                       cw,Hml,K,aml,L,hc,
                       a,b,c,d,
                       Ka,Kb,Kc,Kd)

    