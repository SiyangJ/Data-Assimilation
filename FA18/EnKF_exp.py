import EnKF_func
import AuxFuncs
import ML_pipeline
import pandas
from tensorflow import keras

if __name__=="__main__":
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_dir",type=str)
    parser.add_argument("data_dir",type=str)
    
    parser.add_argument("--RSEED",type=int,default=522)
    parser.add_argument("--sigmaobs",type=float,default=0.9)
    parser.add_argument("--nobs",type=int,default=200)
    parser.add_argument("--HML",type=bool,default=False)
    parser.add_argument("--model_dir",type=str,default=None)
    
    parser.add_argument("--infl_nlin",type=float,default=1.0)
    parser.add_argument("--sigmainit",type=float,default=1.3)
    parser.add_argument("--nens",type=int,default=100)
    parser.add_argument("--ferr",type=float,default=1.0)
    parser.add_argument("--deltaobs",type=float,default=0.1)
    

    
    args = parser.parse_args()
    
    H = None
    HTYPE = 0
    
    if args.HML:
        HTYPE = None
        m = keras.models.load_model(args.model_dir)
        H = ML_pipeline.MakeHML(m)
    
    data_dir = args.data_dir+'da{:.3f}.npz'.format(time.time())
    
    f_all,a_all,f_10, a_10,f_30, a_30 =EnKF_func.EnKF(
                               RSEED=args.RSEED,
                               sigmaobs=args.sigmaobs,
                               nobs=args.nobs,
                               SAVEDATA=3,
                               data_dir=data_dir,
                               HTYPE=HTYPE,
                               H=H,
                               infl_nlin=args.infl_nlin,
                               sigmainit=args.sigmainit,
                               nens=args.nens,
                               ferr=args.ferr,
                               deltaobs=args.deltaobs)
    
    import pandas as pd
    df = pd.read_csv(args.csv_dir)
    df.loc[df.shape[0]] = ([data_dir, model_dir,
                            args.RSEED,     args.sigmaobs, args.nobs, args.deltaobs,
                            args.sigmainit, args.nens,     args.ferr, args.infl_nlin,
                            f_all, a_all, f_10, a_10, f_30, a_30])
    df.to_csv(args.csv_dir,index=False)
    