import ML_pipeline

if __name__=="__main__":
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("RSEED",type=int)
    parser.add_argument("sigmaobs",type=float)
    parser.add_argument("pickle_dir",type=str)
    
    args = parser.parse_args()
    
    mae,acc=ML_pipeline.ML_exp(RSEED=args.RSEED,sigmaobs=args.sigmaobs,nobs=2000)
    
    import pickle
    import numpy as np
    f = open(args.pickle_dir,'rb+') 
    data_arr = pickle.load(f)
    f.close()
    f = open(args.pickle_dir,'wb+') 
    print(data_arr)
    to_append = [[args.RSEED,args.sigmaobs,mae,acc]]
    print(to_append)
    np.append(data_arr,to_append,axis=0)
    print(data_arr)
    pickle.dump(data_arr,f)
    f.close()
    