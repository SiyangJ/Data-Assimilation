import numpy as np
from config import CFP

import ObservationOperator
import Model

def DataGeneration():
    
    random_seed = CFP['DataGeneration'].getint('random_seed')
    np.random.seed(random_seed)
    ndim = CFP['TrueModel'].getint('dimension')
    # True model
    
    M = getattr(Model,CFP['TrueModel']['model_equation'])
    
    resume_from_last = CFP['DataGeneration'].getboolean('resume_from_last')
    if resume_from_last:
        last_data = np.load(CFP['DataGeneration']['last_dir'])
    
    if resume_from_last and CFP['DataGeneration'].getboolean('tinit'):
        tinit = last_data['tinit']
    else:
        tinit = 0
    if resume_from_last and CFP['DataGeneration'].getboolean('xinit'):
        xinit = last_data['xinit']
    else:
        xinit = np.random.rand(ndim)
        if CFP['DataGeneration'].getboolean('from_attractor'):
            for i in range(100):
                xinit = M(xinit,i,i+1)
            tinit = 100
    
    if resume_from_last and truestate:
        truestate = last_data['truestate']
    else:
        num_data = CFP['DataGeneration'].getint('num')
        delta_data = CFP['DataGeneration'].getfloat('delta')
        truestate = np.zeros([num_data,ndim])
        truestate[0,:] = xinit
        ttemp = tinit
        for i in range(1,num_data):
            truestate[i,:] = M(truestate[i-1,:],ttemp,ttemp + delta_data)
            ttemp += delta_data
    
    np.savez(CFP['DataGeneration']['save_dir'],
             xinit = xinit,
             tinit = tinit,
             random_seed=random_seed,
             truestate=truestate)
    
def ObservationGeneration():
    H = getattr(ObservationOperator,CFP['ObservationGeneration']['observation_operator'])
    if STAGE < 3 or trueobs is None:
        ## generate observations
        trueobs = H(truestate)
    
    robs = sigmaobs*sigmaobs
    Robsmat = np.identity(dobs)*robs
    
    if yobs is None:
        yobs = trueobs + np.random.multivariate_normal(np.zeros(dobs), Robsmat,nobs+1).T