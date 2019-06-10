import numpy as np
import config
import Model
import ObservationOperator

CFP = config.CFP

class EnKF():
    
    def __init__(self):
        # True model
        self.TM = getattr(Model,CFP['TrueModel']['model_equation'])
        # Forcast model
        self.TF = getattr(Model,CFP['RunModel']['model_equation'])
        self.Hm = getattr(ObservationOperator,CFP['DataAssimilation']['observation_operator'])
        
        self.linear = CFP['DataAssimilation'].getboolean('linear')
        self.random_seed = CFP['DataAssimilation'].getfloat('random_seed')
        
        self.number_observation = CFP['DataGeneration'].getint('num')
        self.dimension_state = CFP['RunModel'].getint('dimension')
        self.dimension_observation = CFP['ObservationGeneration'].getint('dimension_observation')
        self.number_ensemble = CFP['DataAssimilation'].getint('number_ensemble')
        
        data = np.load(CFP['DataAssimilation']['data_source'])
        obs  = np.load(CFP['DataAssimilation']['observation_source'])
        
        self.truestate = data['truestate']
        self.tinit = data['tinit']
        self.tdelta = data['tdelta']
        self.tobs = np.range(self.number_observation) * self.tdelta + tinit
        
        self.trueobs = obs['true_observation']
        self.yobs = obs['output_observation']
        
        self.sigma_init = CFP['DataAssimilation'].getfloat('sigma_init')
        
        self.inflation = CFP['DataAssimilation'].getfloat('inflation')
        
        self.value_bound_upper = config.GetArray('RunModel','value_bound_upper')
        self.value_bound_lower = config.GetArray('RunModel','value_bound_lower')

    def EnKF(self):
        
        nobs = self.number_observation
        dobs = self.dimension_observation
        ndim = self.dimension_state
        nens = self.number_ensemble
        
        ## Containers
        xfm = np.zeros([nobs, ndim])
        xam = np.zeros([nobs, ndim])
        xfcov = np.zeros([nobs, ndim, ndim])
        xacov = np.zeros([nobs, ndim, ndim])

        yfm = np.zeros([nobs, dobs])
        yfcov = np.zeros([nobs, dobs, dobs])
        
        xinit = self.truestate[0,:]
        sigmainit = self.sigma_init
        
        ## set the initial mean and covariance:
        
        sigma = xinit * np.clip(sigmainit * np.random.randn(*xinit.shape),-0.5,1)
        
        xfm[0,:] = xam[0,:] = xinit + sigma
        xfcov[0,:,:] = xacov[0,:,:] = np.identity(ndim) * (sigma ** 2)
        
        ## Generate initial ensemble
        xens = np.random.multivariate_normal(xam[0,:], xacov[0,:,:], nens)
        
        xens = np.clip(xens,self.value_bound_lower,self.value_bound_upper)

        ## Iterate for each observation
        for ii in np.arange(1,nobs):

            ## Forecast
            for jj in np.arange(nens):
                xens[jj,:] = TF(xens[jj,:], self.tobs[ii-1], self.tobs[ii])

            ## Forecast mean and cov
            xfm[ii,:] = np.mean(xens,axis=0)
            xfcov[ii,:,:] = np.cov(xens.T)

            ## Observations generated from forecast
            ## Estimating the covariances with sample covariance
            yens = self.Hm(xens)
            yfm[ii,:] = np.mean(yens,axis=0)
            yfcov[ii,:,:] = np.cov(yens.T)

            ## Kalman gain:
            if self.linear:
            # Linear EnKF
                pf = xfcov[ii,:,:]
                
                ## TODO
                
                if INFLATION:
                # inflation
                    pf += infl_lin*np.trace(xfcov[ii,:,:])/ndim*np.identity(ndim)

                pfht = pf @ Hmat.T
                kgain = pfht @ np.linalg.pinv( Hmat @ pfht + Robsmat)

            else:
            # Nonlinear EnKF
                Ex = xens - xfm[ii,:]
                Ey = yens - yfm[ii,:]    
                Bxy = Ex.T @ Ey / (nens-1)
                Byy = yfcov[ii,:,:]

                # Inflation
                Byy += self.inflation * Robsmat
                kgain = Bxy @ np.linalg.pinv(Byy)

            ## Update ensemble members
            pertobs = np.random.multivariate_normal(yobs[:,ii],Robsmat,nens)

            for jj in np.arange(nens):
                xens[jj,:] += kgain @ (pertobs[jj,:] - Hm(xens[[jj],:].T).flatten())

            ## Analysis mean and cov
            xam[ii,:] = np.mean(xens,axis=0)
            xacov[ii,:,:] = np.cov(xens.T)


        ## Calculate statistics
        # Absolute errors
        xferr = np.abs(xfm - truestate.T)
        xaerr = np.abs(xam - truestate.T)

        xferr10 = xferr[:,0:10]
        xferr30 = xferr[:,10:]
        xaerr10 = xaerr[:,0:10]
        xaerr30 = xaerr[:,10:]

        xferr10avgx = np.mean(xferr10,axis=1)
        xaerr10avgx = np.mean(xaerr10,axis=1)

        xferr30avgx = np.mean(xferr30,axis=1)
        xaerr30avgx = np.mean(xaerr30,axis=1)

        xferravgx = np.mean(xferr,axis=1)
        xaerravgx = np.mean(xaerr,axis=1)

        avgstp = 30
        xferravgt = np.mean(xferr[avgstp:,:],axis=0)
        xferravgxt = np.mean(xferravgt)

        xaerravgt = np.mean(xaerr[avgstp:,:],axis=0)
        xaerravgxt = np.mean(xaerravgt)

        xferr10avgxt = np.mean(xferr10avgx)
        xaerr10avgxt = np.mean(xaerr10avgx)

        xferr30avgxt = np.mean(xferr30avgx)
        xaerr30avgxt = np.mean(xaerr30avgx)

        np.savez(data_dir,
                    ## Paramters
                    # Notebook usage flags
                    DESC=DESC,
                    RSEED=RSEED,
                    HTYPE=HTYPE,
                    HERROR=HERROR,
                    LINEAR=LINEAR,
                    FORECAST_ERROR=FORECAST_ERROR,
                    INFLATION=INFLATION,
                    # Model parameters
                    ndim=ndim,
                    pars=pars,
                    # Observation paramters
                    nobs=nobs,
                    deltaobs=deltaobs,
                    dobs=dobs,
                    Hmat=Hmat,
                    sigmaobs=sigmaobs,

                    ## DA paramters
                    infl_lin=infl_lin,
                    infl_nlin=infl_nlin,
                    sigmainit=sigmainit,
                    nens=nens,
                    ferr=ferr,

                    ## Truth
                    truestate=truestate,
                    trueobs=trueobs,
                    yobs=yobs,
                    ## Running Data
                    xfm=xfm,xfcov=xfcov,
                    xam=xam,xacov=xacov,
                    yfm=yfm,yfcov=yfcov,
                    ## Statistics
                    xferravgxt  =xferravgxt,  xaerravgxt  =xaerravgxt,
                    xferr10avgxt=xferr10avgxt,xaerr10avgxt=xaerr10avgxt,
                    xferr30avgxt=xferr30avgxt,xaerr30avgxt=xaerr30avgxt)

        return (xferravgxt,   xaerravgxt,
                xferr10avgxt, xaerr10avgxt,
                xferr30avgxt, xaerr30avgxt)