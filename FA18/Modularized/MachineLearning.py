import numpy as np
import scipy.misc
import h5py
import pandas as pd
import tensorflow as tf
import time
import random
import numpy as np
import os
import sys
from copy import deepcopy
from config import CFP

def prepare_dirs(delete_train_dir=False):
    # Create checkpoint dir (do not delete anything)
    if not tf.gfile.Exists(CFP['MachineLearning']['save_dir']):
        tf.gfile.MakeDirs(CFP['MachineLearning']['save_dir'])
    
    # Cleanup train dir
    if delete_train_dir:
        if tf.gfile.Exists(CFP['MachineLearning']['save_dir']):
            tf.gfile.DeleteRecursively(CFP['MachineLearning']['save_dir'])
        tf.gfile.MakeDirs(CFP['MachineLearning']['save_dir'])

def setup_tensorflow_test():
    
    tfconfig = tf.ConfigProto()
    sess = tf.Session(config=tfconfig)

    # Initialize rng with a deterministic seed
    with sess.graph.as_default():
        tf.set_random_seed(CFP['MachineLearning'].getint('random_seed'))
        
    random.seed(CFP['MachineLearning'].getint('random_seed'))
    np.random.seed(CFP['MachineLearning'].getint('random_seed'))

    tf.gfile.MkDir('%s/test_log' % (CFP['MachineLearning']['save_dir'],))
    summary_writer = tf.summary.FileWriter('%s/test_log' % (CFP['MachineLearning']['save_dir'],), sess.graph)
    return sess, summary_writer

def setup_tensorflow_train():
    
    tfconfig = tf.ConfigProto()
    sess = tf.Session(config=tfconfig)

    # Initialize rng with a deterministic seed
    with sess.graph.as_default():
        tf.set_random_seed(CFP['MachineLearning'].getint('random_seed'))
        
    random.seed(CFP['MachineLearning'].getint('random_seed'))
    np.random.seed(CFP['MachineLearning'].getint('random_seed'))

    tf.gfile.MkDir('%s/training_log' % (CFP['MachineLearning']['save_dir'],))
    tf.gfile.MkDir('%s/validation_log' % (CFP['MachineLearning']['save_dir'],))
    summary_writer = tf.summary.FileWriter('%s/training_log' % (CFP['MachineLearning']['save_dir'],), sess.graph)
    val_sum_writer = tf.summary.FileWriter('%s/validation_log' % (CFP['MachineLearning']['save_dir'],), sess.graph)

    return sess, summary_writer, val_sum_writer

def save_checkpoint(train_data,batch):
    td = train_data
    saver = tf.train.Saver()
    model_path = os.path.join(CFP['MachineLearning']['save_dir'],'checkpoint_'+str(batch))
    
    save_path = saver.save(td.sess, model_path)
    print ('Model saved in file: %s' % saver.last_checkpoints)

class TrainData(object):
    def __init__(self, dictionary):
        self.__dict__.update(dictionary)
        
class TestData(object):
    def __init__(self, dictionary):
        self.__dict__.update(dictionary)
        
class DataGenerator:
    def __init__(self,for_evaluation = None):        
        if for_evaluation is None:
            for_evaluation = sys.argv[1] in ['inference','evaluation']
        
        self.for_evaluation = for_evaluation       
        self._X,self._Y = DataGenerator._OrganizeData(self.for_evaluation)
        
        self.normalization = CFP['MachineLearning'].getboolean('normalization',False)
        self.sep = CFP['MLTraining'].getfloat('validation_proportion')
        
        _res = DataGenerator._PrepareData(self._X,self._Y,
                            normalization = self.normalization,
                            sep = self.sep)
        self._train_X, self._train_Y, self._test_X, self._test_Y = _res[:4]

        if self.normalization:
            self._Xm, self._Xs, self._Ym, self._Ys = _res[4:8]
            np.savez(os.path.join(CFP['MachineLearning']['save_dir'],'normalization_params.npz'),
                     Xm=self._Xm, Xs=self._Xs, Ym=self._Ym, Ys=self._Ys)
    
    def get_training_and_testing_generators(self):        
        batch_size = CFP['MLTraining'].getint('batch_size')
        training_generator = DataGenerator.data_generator(self._train_X, self._train_Y, batch_size=batch_size)
        validation_generator = DataGenerator.data_generator(self._test_X, self._test_Y, batch_size=batch_size)
        return training_generator, validation_generator
    
    def _OrganizeData(for_evaluation = False):
        X = np.empty([0,CFP['MachineLearning'].getint('input_dim')])
        Y = np.empty([0,CFP['MachineLearning'].getint('output_dim')])

        CFP_sec = CFP['MLEvaluation'] if for_evaluation else CFP['MLTraining']
        data_dir = CFP_sec['data_source']
        obs_dir = CFP_sec['observation_source']

        for f in os.listdir(obs_dir):
            f_full = os.path.join(obs_dir,f)
            if os.path.isfile(f_full) and (f.find('.npz')>=0 or f.find('.npy')>=0):
                try:
                    Yfile = np.load(f_full)
                    Y = np.append(Y,Yfile['output_observation'],axis=0)
                    Xfile = np.load(str(Yfile['data_dir']))
                    X = np.append(X,Xfile['truestate'],axis=0)
                except:
                    continue
        return X,Y

    def _PrepareData(_X,_Y,normalization=False,sep=0.8):    
        L = _X.shape[0]
        
        train_size=int(L * (1-sep))
        index = np.arange(L)
        np.random.shuffle(index)
        train_index = index[:train_size]
        test_index = index[train_size:]

        X = _X[train_index,:]
        Y = _Y[train_index,:]
        Xt = _X[test_index,:]
        Yt = _Y[test_index,:]
        if not normalization:
            return X,Y,Xt,Yt
        else:
            Xm = X.mean(axis=0)
            Ym = Y.mean(axis=0)
            Xs = X.std(axis=0)
            Ys = Y.std(axis=0)
            return (X-Xm)/Xs,(Y-Ym)/Ys,(Xt-Xm)/Xs,(Yt-Ym)/Ys,Xm,Xs,Ym,Ys

    def data_generator(X,Y,batch_size=0):
        L = X.shape[0]
        i = 0
        ## TODO
        ## batch method
        while True:
            if batch_size == 0:
                yield X,Y
                continue
            if i == L:
                yield None,None
                i = 0
                continue
            yield X[[i,],:],Y[[i,],:]
            i += 1      

class MachineLearning:
    def __init__(self,for_evaluation = None):
        if for_evaluation is None:
            for_evaluation = sys.argv[1] in ['inference','evaluation']
        self.for_evaluation = for_evaluation       
        self.dataGen = DataGenerator(for_evaluation)

    def cal_loss(pred,truth):
        loss=tf.reduce_mean(tf.square(pred-truth))
        return loss

    def create_model(train_phase=True):
        x_size = CFP['MachineLearning'].getint('input_dim')
        y_size = CFP['MachineLearning'].getint('output_dim')

        input_shape = [None,x_size]    
        output_shape = [None,y_size]

        old_vars = tf.global_variables()

        X_variable = tf.placeholder(tf.float32, shape=input_shape, name='X_placeholder')
        Y_variable = tf.placeholder(tf.float32, shape=output_shape, name='Y_placeholder')

        keep_prob = CFP['MLTraining']['keep_prob']

        layers_string = CFP['MachineLearning']['layers']    
        layers_strings = layers_string.split(',')    
        layers = [int(s) for s in layers_strings]

        cur_layer = X_variable
        for layer_num in layers:
            cur_layer = tf.layers.dense(cur_layer,layer_num,
                                        activation = tf.nn.relu,
                                        bias_initializer=tf.constant_initializer(0.01))
            cur_layer = tf.layers.dropout(cur_layer,rate = keep_prob)

        ## Default is linear
        pred = tf.layers.dense(cur_layer,y_size,bias_initializer=tf.constant_initializer(0.01))

        new_vars = tf.global_variables()
        gene_vars = list(set(new_vars)-set(old_vars))

        loss = MachineLearning.cal_loss(pred,Y_variable)
        tf.summary.scalar('loss',loss)

        final_loss = loss
        # l2 loss
        # for each varibale name, like:  conv2a/biases:0  and aux2_pred/weights:0
        if CFP['MLTraining'].getboolean('regularize'):
            for _var in gene_vars:
                # to use L2 loss, all the variables must be with the name of "weights"
                if _var.name.find('kernel') > -1 or _var.name.find('weight') > -1: # to exclude 'bias' term
                    final_loss = final_loss + CFP['MLTraining'].getfloat('reg_lambda') * tf.nn.l2_loss(_var)

        return (X_variable, Y_variable,
                pred, loss, final_loss,
                gene_vars)

    def create_optimizers(train_loss):
        learning_rate  = tf.placeholder(dtype=tf.float32, name='learning_rate')
        train_opti = tf.train.AdamOptimizer(learning_rate)
        global_step    = tf.Variable(0, dtype=tf.int64, trainable=False, name='global_step')
        # Disable var_list, train all variables
        # var_list = _get_var_list()
        train_minimize = train_opti.minimize(train_loss, name='loss_minimize', global_step=global_step)
        return train_minimize, learning_rate, global_step

    def train_model(self,train_data):
        td = train_data
        summaries = tf.summary.merge_all()

        if CFP['MLTraining'].getboolean('restore_from_last'):
            saver = tf.train.Saver()
            model_path = tf.train.latest_checkpoint(CFP['MLTraining']['last_trained_checkpoint'])
            print('Machine Learning Training: restore last checkpoint from:%s' % model_path)
            saver.restore(td.sess, model_path)
        else:
            init_op = tf.global_variables_initializer()
            print('Machine Learning Training: global variable initialization...')
            td.sess.run(init_op)

        lrval       = CFP['MLTraining'].getfloat('initial_learning_rate')
        start_time  = time.time()
        done  = False
        epoch = 0    
        training_generator, testing_generator = self.dataGen.get_training_and_testing_generators()

        decay_after = CFP['MLTraining'].getint('decay_after')
        decay_ratio = CFP['MLTraining'].getfloat('decay_ratio')
        max_epoch = CFP['MLTraining'].getint('epochs')
        save_every_n = CFP['MLTraining'].getint('save_every_n')
        summary_every_n = CFP['MLTraining'].getint('summary_every_n')

        while not done:
            epoch += 1        
            train_X, train_Y = next(training_generator)
            feed_dict = {td.X_variable : train_X, 
                         td.Y_variable : train_Y, 
                         td.learning_rate : lrval}

            cur_time = time.ctime()

            if epoch % summary_every_n == 0:
                ops = [td.train_minimize, td.loss, summaries,] 
                [_, loss, summary_vis] = td.sess.run(ops, feed_dict=feed_dict)

                print('[%25s], epoch [%4d], lr[%1.8f] ,loss[%3.10f]'
                                % (cur_time, epoch, lrval, loss) )

                # Update learning rate
                if epoch % decay_after == 0:
                    lrval *= decay_ratio

                td.summary_writer.add_summary(summary_vis, epoch)

                val_X, val_Y = next(testing_generator)            
                val_feed_dict = {td.X_variable : val_X, 
                                 td.Y_variable : val_Y}

                val_ops = [td.loss, summaries,] 
                [val_loss, val_summary] = td.sess.run(val_ops, feed_dict=val_feed_dict)

                print('[%25s], validation: iter [%4d], loss[%3.10f]'
                                % (cur_time, epoch, val_loss) )

                td.val_sum_writer.add_summary(val_summary, epoch)

            else:
                ops = [td.train_minimize, td.loss] 
                [_, loss] = td.sess.run(ops, feed_dict=feed_dict)

            if epoch % save_every_n == 0:
                save_checkpoint(td, epoch)

            if epoch  >= max_epoch:
                done = True
        save_checkpoint(td, epoch)
        print('Finished training!')

    def train(self):
        prepare_dirs(delete_train_dir=False)
        sess, summary_writer, val_sum_writer = setup_tensorflow_train()

        (X_variable, Y_variable,
         pred, loss, final_loss,
         gene_vars) = MachineLearning.create_model(train_phase=True)

        train_minimize, learning_rate, global_step = MachineLearning.create_optimizers(final_loss)

        train_data = TrainData(locals())
        MachineLearning.train_model(train_data)

    def evaluate(self,td,X=None,Y=None):
        start_time  = time.time()
        if X is None:
            X = self.dataGen._X
        if Y is None:
            Y = self.dataGen._Y
        feed_dict = { td.X_variable : X, 
                     td.Y_variable : Y}
        ops = [td.pred,td.loss]
        [pred,loss] = td.sess.run(ops, feed_dict=feed_dict)
        elapsed = int(time.time() - start_time)
        print('evaluate complete, cost [%3d] seconds' % (elapsed))

        return pred,loss  

    def test(self):
        if self.dataGen.normalization:
            norm_params = np.load(os.path.join(CFP['MachineLearning']['save_dir'],'normalization_params.npz'))
            X = (self.dataGen._X - norm_params['Xm']) / norm_params['Xs']
            Y = (self.dataGen._Y - norm_params['Ym']) / norm_params['Ys']
        else:
            X = self.dataGen._X
            Y = self.dataGen._Y

        prepare_dirs(delete_train_dir=False)
        sess, summary_writer = setup_tensorflow_test()

        (X_variable, Y_variable,
         pred, loss, final_loss,
         gene_vars) = MachineLearning.create_model(train_phase=True)

        saver = tf.train.Saver()
        model_path = tf.train.latest_checkpoint(CFP['MachineLearning']['save_dir'])
        print('saver restore from:%s' % model_path)
        saver.restore(sess, model_path)

        test_data = TestData(locals())
        pred,loss = MachineLearning.evaluate(test_data,X=X,Y=Y)

        if self.dataGen.normalization:
            X = X * norm_params['Xs'] + norm_params['Xm']
            Y = Y * norm_params['Ys'] + norm_params['Ym']
            pred = pred * norm_params['Ys'] + norm_params['Ym']

        mse_loss = ((pred-Y)**2).mean()
        evals = {'MSE':mse_loss,}
        if CFP['MLEvaluation'].getboolean('save_predict'):
            np.savez(CFP['MLEvaluation']['save_dir'],X=X,Y=Y,predict=pred,evaluation=evals)
        return pred,mse_loss

def main(argv=None):
    train()

if __name__ == '__main__':
    tf.app.run()