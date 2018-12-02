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

class TrainData(object):
    def __init__(self, dictionary):
        self.__dict__.update(dictionary)
        
class TestData(object):
    def __init__(self, dictionary):
        self.__dict__.update(dictionary)

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
    
    loss = cal_loss(pred,Y_variable)
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
    global_step    = tf.Variable(0, dtype=tf.int64,   trainable=False, name='global_step')
    # Disable var_list, train all variables
    # var_list = _get_var_list()
    train_minimize = train_opti.minimize(train_loss, name='loss_minimize', global_step=global_step)
    return train_minimize, learning_rate, global_step

def _OrganizeData():
    
    X = np.empty([0,CFP['MachineLearning'].getint('input_dim')])
    Y = np.empty([0,CFP['MachineLearning'].getint('output_dim')])
    
    for_evaluation = sys.argv[1] in ['inference','evaluation']
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
    #X = np.array(X).astype('float32')
    #Y = np.array(Y).astype('float32')
    return X,Y

_X,_Y = _OrganizeData()

def _PrepareData():    
    L = _X.shape[0]
    sep = CFP['MLTraining'].getfloat('validation_proportion')
    
    train_size=int(L * (1-sep))
    index = np.arange(L)
    np.random.shuffle(index)
    train_index = index[:train_size]
    test_index = index[train_size:]
    
    return _X[train_index,:],_Y[train_index,:],_X[test_index,:],_Y[test_index,:]

_train_X,_train_Y,_test_X,_test_Y = _PrepareData()

def data_generator(istrain=True,batch_size=0):
    X = _train_X if istrain else _test_X
    Y = _train_Y if istrain else _test_Y
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

def get_training_and_testing_generators():
    batch_size = CFP['MLTraining'].getint('batch_size')
    training_generator = data_generator(istrain=True, batch_size=batch_size)
    validation_generator = data_generator(istrain=False, batch_size=batch_size)
    
    return training_generator, validation_generator

def _save_checkpoint(train_data,batch):
    td = train_data
    saver = tf.train.Saver()
    model_path = os.path.join(CFP['MachineLearning']['save_dir'],'checkpoint_'+str(batch))
    
    save_path = saver.save(td.sess, model_path)
    print ('Model saved in file: %s' % saver.last_checkpoints)

def train_model(train_data):
    td = train_data
    summaries = tf.summary.merge_all()

    # if the sess is restored from last checkpoint, do not need to 
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
    training_generator, testing_generator = get_training_and_testing_generators()

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
            _save_checkpoint(td, epoch)

        if epoch  >= max_epoch:
            done = True

    _save_checkpoint(td, epoch)
    print('Finished training!')

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

def prepare_dirs(delete_train_dir=False):
    # Create checkpoint dir (do not delete anything)
    if not tf.gfile.Exists(CFP['MachineLearning']['save_dir']):
        tf.gfile.MakeDirs(CFP['MachineLearning']['save_dir'])
    
    # Cleanup train dir
    if delete_train_dir:
        if tf.gfile.Exists(CFP['MachineLearning']['save_dir']):
            tf.gfile.DeleteRecursively(CFP['MachineLearning']['save_dir'])
        tf.gfile.MakeDirs(CFP['MachineLearning']['save_dir'])

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

def train():
    prepare_dirs(delete_train_dir=False)
    sess, summary_writer, val_sum_writer = setup_tensorflow_train()

    (X_variable, Y_variable,
     pred, loss, final_loss,
     gene_vars) = create_model(train_phase=True)

    train_minimize, learning_rate, global_step = create_optimizers(final_loss)

    train_data = TrainData(locals())
    train_model(train_data)
    
## A very crude version.
## Just for experiment.
def predict(td):
    start_time  = time.time()
    feed_dict = { td.X_variable : _X, 
                 td.Y_variable : _Y}
    ops = [td.pred,td.loss]
    [pred,loss] = td.sess.run(ops, feed_dict=feed_dict)
    elapsed = int(time.time() - start_time)
    print('Predict complete, cost [%3d] seconds' % (elapsed))

    return pred,loss  

def test():
    prepare_dirs(delete_train_dir=False)
    sess, summary_writer = setup_tensorflow_test()

    (X_variable, Y_variable,
     pred, loss, final_loss,
     gene_vars) = create_model(train_phase=True)

    saver = tf.train.Saver()
    model_path = tf.train.latest_checkpoint(CFP['MachineLearning']['save_dir'])
    print('saver restore from:%s' % model_path)
    saver.restore(sess, model_path)
    
    test_data = TestData(locals())
    pred,loss = predict(test_data)
    evals = {'MSE':loss,}
    if CFP['MLEvaluation'].getboolean('save_predict'):
        np.savez(CFP['MLEvaluation']['save_dir'],X=_X,Y=_Y,predict=pred,evaluation=evals)
    return pred,loss

def main(argv=None):
    train()

if __name__ == '__main__':
    tf.app.run()