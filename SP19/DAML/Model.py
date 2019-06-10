import numpy as np
import math
import config

def EW09default(x0,t0,t1):
    
    if isinstance(x0,np.ndarray):
        x0 = x0.tolist()
    
    import matlab as ml

    eng = config.MATLAB_ENGINE
    Fc = config.CFP['TrueModel'].getfloat('Fc')
    return np.array(eng.Mdefault(eng.transpose(ml.double(x0)),float(t0),float(t1),float(Fc)))[0]

def EW09(x0,t0,t1):
    if isinstance(x0,np.ndarray):
        x0 = x0.tolist()
    
    import matlab as ml

    eng = config.MATLAB_ENGINE
    
    M = config.CFP['TrueModel']    
    Fc = M.getfloat('Fc')
    Fs = M.getfloat('Fs')
    Fo = M.getfloat('Fo')
    Ft = M.getfloat('Ft')
    Fb = M.getfloat('Fb')
    cw = M.getfloat('cw')
    Hml = M.getfloat('Hml')
    K = M.getfloat('K')
    aml = M.getfloat('aml')
    L = M.getfloat('L')
    hc = M.getfloat('hc')
    a = M.getfloat('a')
    b = M.getfloat('b')
    c = M.getfloat('c')
    d = M.getfloat('d')
    Ka = M.getfloat('Ka')
    Kb = M.getfloat('Kb')
    Kc = M.getfloat('Kc')
    Kd = M.getfloat('Kd')
    
    return np.array(eng.M(eng.transpose(ml.double(x0)),float(t0),float(t1),float(Fc),
                       Fs,Fo,Ft,Fb,
                       cw,Hml,K,aml,L,hc,
                       a,b,c,d,
                       Ka,Kb,Kc,Kd))[0]