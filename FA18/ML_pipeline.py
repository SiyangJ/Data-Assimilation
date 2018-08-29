
# coding: utf-8

# ## Summary:
# **Name:** Functionalize ML pipeline  
# **Author:** Siyang Jing  
# **Organization:** UNC-CH  
# **License:** WTFPL  
# 
# **Reference:**
# 1. Tensorflow and Keras
# 1. The author's other codes
# 1. Relevant numerous papers
#   
# **Description:**  
# This file prepares a function ML pipeline based on ML_pipeline.ipynb. This code sets up the numerical experiment with only machine learning part. Input parameters of model, observation operator, ML architecture, range of value, and receive statistics.
# _Input_:
# * Parameters: (To be implemented)
#   * ndim: dimension of state space default is 40,
#   * pars: parameters for model ode, default value is (8.0) for Lorenz 96,
#   * nobs: number of observations, dafault is 50,
#   * deltaobs=0.1,
#   * dobs=None,
#   * Hmat=None,
#   * sigmaobs=0.9,
#   * infl_lin=0.5,
#   * infl_nlin=1,
#   * sigmainit=1.3,
#   * nens=100,
#   * ferr=1.0
# * DA Flags: (To be implemented)
#   * HTYPE: type of observation operator
#     * None: everything is specified in parameters, **TODO**
#     * 0: ObsOp_40_20
#     * 1: ObsOp_40_20 with Inv_20_10
#   * LINEAR: whether to use linearized KF or not, specific to $H$
#   * HERROR: whether H and Hm are the same
# * Usage Flags: (To be implemented)
#   * PLOTTING
#   * DEBUG
#   * SAVEDATA: indicator for which data to save  
#     * 0: Don't save anything
#     * 1: Save everything including,  
#       1. Flags and Parameters  
#       1. xfm, xam, xfcov, xacov  
#       1. yfm, yfcov  
#       1. statistics  
#     * 2: Flags, Paramters, Statistics (No running data)
#     * 3: haha
#   
# _Output_:  
# Forecast/analysis absolute error averaged over variables and time (from the 30th observation step as the starting point after DA stabilizes)
# 1. xferravgxt: forecast absolute error averaged over all variables and time
# 1. xaerravgxt: analysis absolute error averaged over all variables and time
# 1. xferr10avgxt: forecast absolute error averaged over first 10 variables and time
# 1. xaerr10avgxt: analysis absolute error averaged over first 10 variables and time
# 1. xferr30avgxt: forecast absolute error averaged over last 30 variables and time
# 1. xaerr30avgxt: analysis absolute error averaged over last 30 variables and time
# 
# _Saved Files_:   
# Variables saved in a .npz file with their orginal names with numpy.savez method.  
# 
# _Plots_:
# 1. Long true trajectory
# 1. Run results
# 1. Statistics
# 
# **Requirements:**
# 1. Relevant Python modules
# 1. AuxFuncs, which defines the following:  
#   1. Observation operator  
#   1. Stupid inverse function  
#   1. Lorenz96 as model

# ## TODO List
# 1. Needs a more sophisticated build_model function
#   * Input:
#     1. Layers
#     1. Optimizer
#   * Architecture selection
# 1. 

from EnKF_func import *
from AuxFuncs import *
from ML_func import *
from DataGenerator import *

TRAIN_SPLIT = 0.8
NORMALIZE = False
LOSS = 'mse'
EPOCHS = 1000

def build_model(loss=LOSS,layers=[20,20]):
    model = keras.Sequential()
    model.add(keras.layers.Dense(layers[0], 
                                 activation=tf.nn.relu, 
                                 input_shape=(40,)))
    for num in layers[1:]:
        model.add(keras.layers.Dense(num, activation=tf.nn.relu))
    model.add(keras.layers.Dense(20))

    model.compile(loss=loss,
                optimizer='adam',
                metrics=['mae','acc'])
    return model

def PrepareData(X,Y,
                train_split=TRAIN_SPLIT,
                normalize=NORMALIZE):
    
    data_size = np.shape(X)[0]

    # Shuffle the input
    order = np.argsort(np.random.random(data_size))

    train_size = round(data_size * train_split)
    test_size = data_size - train_size

    train_data = X[order[0:train_size],:]
    train_labels = Y[order[0:train_size],:]

    test_data = X[order[train_size:],:]
    test_labels = Y[order[train_size:],:]

    if normalize:
        ## Normalize the data
        ## Doesn't seem particularly useful
        # Test data is *not* used when calculating the mean and std.
        '''
        mean = train_data.mean(axis=0)
        std = train_data.std(axis=0)
        train_data = (train_data - mean) / std
        test_data = (test_data - mean) / std
        '''
    
    return (train_data,train_labels,
            test_data,test_labels)

## Used for model trained with normalized data
'''
def MakeHML(model,std,mean):
    def HML(x):
        x = x.T
        x = (x-mean)/std
        return model.predict(x).T
    return HML
'''
def MakeHML(model):
    def HML(x):
        return model.predict(x.T).T
    return HML


def ML_exp(
    RSEED=215,
    sigmaobs=0.9,
    nobs=1000,
    SAVEDATA=False,
    DATADIR=None,
    SAVEMODEL = False,
    MODELDIR = None,
    EPOCHS=EPOCHS,
    train_split=TRAIN_SPLIT,
    normalize=NORMALIZE,
    loss=LOSS,
    layers=[20,20]
    ):
    X,Y,Y_noise = DataGen(nobs=nobs,RSEED=RSEED,sigmaobs=sigmaobs)

    # Experiment with noised data

    XX = np.transpose(X)
    YY = np.transpose(Y_noise)

    train_data,train_labels,test_data,test_labels = PrepareData(XX,YY,
                                                                train_split=train_split,
                                                                normalize=normalize)

    model1 = build_model(layers=layers,loss=loss)

    m,b,c=ML(train_data,train_labels,test_data,test_labels,
             model=model1,
             EPOCHS=EPOCHS,
             MLDEBUG=False,
             MLPLOTTING=False)
    
    if SAVEMODEL:
        m.save(MODELDIR)
    
    #print('.',end='')
    
    return b,c


def main(RSEED_range,sigmaobs_range):
    RSEED_num = len(RSEED_range)
    sigmaobs_num = len(sigmaobs_range)

    mae_arr = np.zeros([RSEED_num,sigmaobs_num])
    acc_arr = np.zeros([RSEED_num,sigmaobs_num])

    for i in range(RSEED_num):
        for j in range(sigmaobs_num):
            mae_arr[i,j],acc_arr[i,j] = ML_exp(RSEED_range[i],sigmaobs_range[j])

    mae_sigmaobs = np.mean(mae_arr,axis=0)
    acc_sigmaobs = np.mean(acc_arr,axis=0)
    
    return mae_sigmaobs,acc_sigmaobs

RSEED_range = range(215,217)
sigmaobs_range = [0,0.5,1]

if __name__=="__main__":
    mae_sigmaobs, acc_sigmaobs = main(RSEED_range,sigmaobs_range)
    print(mae_sigmaobs,acc_sigmaobs)

