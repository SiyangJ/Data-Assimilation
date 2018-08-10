
# coding: utf-8

# ## Summary:
# **Name:** Data Generation
# **Author:** Siyang Jing  
# **Organization:** UNC-CH  
# **License:** WTFPL  
# 
# **Reference:**
# 1. Amit Apte's EnKF code in Python
# 2. Christian Sampson's EnKF code in matlab
# 3. The author's other codes, esp. EnKF_func.py
# 4. Relevant numerous papers
#   
# **Description:**  
# This file prepares a function EnKF based on Full_EnKF.ipynb.  
# _Input_:
# * Parameters:
#   * ndim: dimension of state space default is 40,
#   * pars: parameters for model ode, default value is (8.0) for Lorenz 96,
#   * nobs: number of observations, dafault is 50,
#   * deltaobs=0.1,
#   * dobs=None,
#   * sigmaobs=0.9,
# * DA Flags:
#   * HTYPE: type of observation operator
#     * None: everything is specified in parameters, **TODO**
#     * 0: ObsOp_40_20
#     * 1: ObsOp_40_20 with Inv_20_10
# * Usage Flags:
#   * SAVEDATA: indicator for whether to save data or not
#   
# _Output_:  
# Truth and observations
# 1. truestate
# 2. trueobs
# 3. yobs
# 
# _Saved Files_:   
# Variables saved in a .npz file with their orginal names with numpy.savez method.  
# 
# **Requirements:**
# 1. Relevant Python modules
# 2. AuxFuncs, which defines the following:  
#   1. Observation operator  
#   2. Stupid inverse function  
#   3. Lorenz96 as model

# ## TODO list:
# 1. Improve the saving procedure
# 2. Write own RK4, standardize the model forward step

import numpy as np
import scipy as scp
import scipy.integrate

import AuxFuncs

def DataGen(
    #### Paramters
    ## Notebook usage flags
    RSEED = 215,
    HTYPE = 0,
    SAVEDATA = False,
    data_dir = None,
    DESC = 'Data',
    ## Model parameters
    ndim = 40,
    pars = (8.0),
    ## Observation paramters
    nobs=50,
    deltaobs=0.1,
    dobs=None,
    H=None,
    # Suggested value
    # Could also choose from long trajectory
    sigmaobs=0.9):
    
    np.random.seed(RSEED)
    
    # True model
    def TM(xin,tin,pars):
        return AuxFuncs.l96rhs(xin,tin,pars)

    if HTYPE is None and H is None:
        if dobs is None:
            dobs = 10
        Hmat = np.zeros([dobs,ndim])
        for i in np.arange(dobs):
            Hmat[i,i] = 1.0            
        def H(x):
            return Hmat @ x
    
    elif HTYPE==0:
        ## Setting for direct observation        
        dobs = 20
        # True observation operator
        def H(x):
            return AuxFuncs.ObsOp_40_20(x)
        
    elif HTYPE==1:
        ## Setting for stupid inversion case
        dobs = 10
        # True observation operator
        def H(x):
            return AuxFuncs.Inv_20_10(AuxFuncs.ObsOp_40_20(x))      

    ## find an initial condition on the attractor

    xrand = np.random.rand(ndim)
    ttrans = np.linspace(0,100,1000)
    xtrans = scp.integrate.odeint(TM, xrand, ttrans, (pars,))

    t = np.linspace(0,50,100000)
    xattr = xtrans[-1,:]

    ## Select sigmaobs based on lone trajectory
    # tlong = np.linspace(0,100,10000)
    # xlong = scp.integrate.odeint(TM, xattr, tlong, (pars,))
    # sigmaobs = np.abs(np.max(xlong[:,0]) - np.min(xlong[:,0]))/50

    ## generate true trajectory

    tend = nobs * deltaobs
    tobs = np.linspace(0, tend, num=nobs+1)
    ttraj = np.linspace(0, tend, num=nobs*100+1)
    truetraj = scp.integrate.odeint(TM, xattr, ttraj, (pars,))
    truestate = truetraj[::100,:]

    truestate = truestate.T

    ## generate observations

    trueobs = H(truestate)
    robs = sigmaobs*sigmaobs
    Robsmat = np.identity(dobs)*robs
    yobs = trueobs + np.random.multivariate_normal(np.zeros(dobs), Robsmat,nobs+1).T
    
    '''
    SAVEDATA: indicator for whether to save data or not
    ''' 
    if SAVEDATA==True:
        np.savez(data_dir,
                ## Paramters
                # Notebook usage flags
                DESC=DESC,
                RSEED=RSEED,
                HTYPE=HTYPE,
                
                # Model parameters
                ndim=ndim,
                pars=pars,
                # Observation paramters
                nobs=nobs,
                deltaobs=deltaobs,
                dobs=dobs,
                sigmaobs=sigmaobs,
                 
                ## Truth
                truestate=truestate,
                trueobs=trueobs,
                yobs=yobs)
    
    return (truestate,trueobs,yobs)