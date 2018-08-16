import ML_pipeline

if __name__=="__main__":
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pickle_dir",type=str)
    parser.add_argument("--RSEED",type=int,default=215)
    parser.add_argument("--sigmaobs",type=float,default=0.9)
    parser.add_argument("--nobs",type=int,default=2000)
    
    args = parser.parse_args()
    
    mae,acc=ML_pipeline.ML_exp(RSEED=args.RSEED,sigmaobs=args.sigmaobs,nobs=args.nobs)
    
    import pickle
    import numpy as np
    f = open(args.pickle_dir,'rb+') 
    data_arr = pickle.load(f)
    f.close()
    f = open(args.pickle_dir,'wb+') 
    #print(data_arr)
    to_append = [[args.RSEED,args.sigmaobs,args.nobs,mae,acc]]
    #print(to_append)
    data_arr = np.append(data_arr,to_append,axis=0)
    #print(data_arr)
    pickle.dump(data_arr,f)
    f.close()
    