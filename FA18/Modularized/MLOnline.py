import numpy as np
import tensorflow as tf
import MachineLearning as ml
import config
from MachineLearning import MachineLearning

CFP = config.CFP

class MLOnline:
    
    def __init__(self):
        self.normalization = CFP['MachineLearning'].getboolean('normalization',False)
        if self.normalization:
            norm_params = np.load(os.path.join(CFP['MachineLearning']['save_dir'],'normalization_params.npz'))
            self._Xm = norm_params['Xm']
            self._Xs = norm_params['Xs']
            self._Ym = norm_params['Ym']
            self._Ys = norm_params['Ys']
        
        tfconfig = tf.ConfigProto()
        self._sess = tf.Session(config=tfconfig)
        with self._sess.graph.as_default():
            tf.set_random_seed(CFP['MachineLearning'].getint('random_seed'))
        
        (self._X_variable,_,
         self._pred,_,_,
         self._gene_vars) = MachineLearning.create_model(train_phase=True)
        
        saver = tf.train.Saver()
        model_path = tf.train.latest_checkpoint(CFP['MachineLearning']['save_dir'])
        print('ML class restores model from:%s' % model_path)
        saver.restore(self._sess, model_path)

    
    def evaluate(self,X):
        if self.normalization:
            X = (X - self._Xm) / self._Xs

        feed_dict = {self._X_variable : X,}
        ops = [self._pred,]
        [pred,] = self._sess.run(ops, feed_dict=feed_dict)

        if self.normalization:
            pred = pred * self._Ys + self._Ym

        return pred