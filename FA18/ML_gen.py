import ML_pipeline
import time

if __name__=="__main__":
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_dir",type=str)
    parser.add_argument("model_dir",type=str)
    parser.add_argument("--RSEED",type=int,default=215)
    parser.add_argument("--sigmaobs",type=float,default=0.9)
    parser.add_argument("--nobs",type=int,default=2000)
    
    args = parser.parse_args()
    
    model_dir = args.model_dir+'m{:.3f}.h5'.format(time.time())
    
    mae,acc=ML_pipeline.ML_exp(RSEED=args.RSEED,
                               sigmaobs=args.sigmaobs,
                               nobs=args.nobs,
                               SAVEMODEL=True,
                               MODELDIR=model_dir)
    
    import pandas as pd
    df = pd.read_csv(args.csv_dir)
    df.loc[df.shape[0]] = ([model_dir, args.RSEED,args.sigmaobs,args.nobs,mae,acc])
    df.to_csv(args.csv_dir)