import numpy as np
import os.path
import scipy.misc
import h5py
import pandas as pd
import tensorflow as tf
import time
import random
import numpy as np
import os
from copy import deepcopy
from config import CFP

def cal_loss(pred,truth):
    loss=tf.reduce_mean(tf.square(pred-truth))
    return loss

def create_model(train_phase=True):
    x_size = CFP['MachineLearning'].getint('input_dim')
    y_size = CFP['MachineLearning'].getint('output_dim')
    
    input_shape = (1,x_size)    
    output_shape = (1,y_size)
    
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

def _PrepareData():
    data_X = pd.read_hdf(FLAGS.data_path,FLAGS.X_ID).values
    data_Y = pd.read_hdf(FLAGS.data_path,FLAGS.Y_ID).values
    
    X,y = [],[]
    total_batch = int(data_X.shape[0] / FLAGS.NUM_PER_DAY)
    for i in range(total_batch):
        X.append(data_X[i*FLAGS.NUM_PER_DAY:(i+1)*FLAGS.NUM_PER_DAY,:])
        y.append(data_Y[i*FLAGS.NUM_PER_DAY:(i+1)*FLAGS.NUM_PER_DAY,:])
    
    train_size=int(FLAGS.sep*len(X))
    split_index=[1]*train_size
    split_index.extend([0] * (len(X) - train_size))
    np.random.shuffle(split_index)

    #division all_data into train and test data
    train_X,train_y,test_X,test_y=[],[],[],[]
    for i,v in enumerate(split_index):
        if v==0:
            test_X.append(X[i])
            test_y.append(y[i])
        else:
            train_X.append(X[i])
            train_y.append(y[i])
    train_X=np.array(train_X).astype('float32')
    train_y=np.array(train_y).astype('float32')
    test_X=np.array(test_X).astype('float32')
    test_y=np.array(test_y).astype('float32')
    return train_X,train_y,test_X,test_y

_train_X,_train_Y,_test_X,_test_Y = _PrepareData()

def get_training_and_testing_generators():
    batch_size = FLAGS.batch_size
    training_generator = data_generator(istrain=True, batch_size=batch_size)
    validation_generator = data_generator(istrain=False, batch_size=batch_size)
    
    return training_generator, validation_generator

def data_generator(istrain=True,batch_size=1):
    X = _train_X if istrain else _test_X
    Y = _train_Y if istrain else _test_Y
    L = X.shape[0]
    S = FLAGS.NUM_PER_DAY
    d = FLAGS.seq_length
    # assert L>=batch_size, "Can't generate batch!"
    while True:
        
        i = np.random.randint(0,L,batch_size)
        j = np.random.randint(0,S-d+1,batch_size)
        
        ret_X,ret_Y=[],[]
        
        for b in range(batch_size):
            ret_X.append(X[i[b]][j[b]:j[b]+d])
            ret_Y.append(Y[i[b]][j[b]:j[b]+d])
            
        ret_X=np.array(ret_X).astype('float32')
        ret_Y=np.array(ret_Y).astype('float32')
        yield ret_X,ret_Y

## A very crude version.
## Just for experiment.
def predict(td):
    preds = []
    X_test = []
    Y_test = []
    start_time  = time.time()
    _, test_generator = get_training_and_testing_generators()
    ## TODO
    ## Just for the purpose of testing.
    patch_num = 20
    print('>> begin predicting for each patch')
    for _i in  range(patch_num):
        _test_X, _test_Y = next(test_generator)
        feed_dict = { td.X_variable : _test_X, 
                     td.Y_variable : _test_Y}
        ops = [td.pred,]
        [pred,] = td.sess.run(ops, feed_dict=feed_dict)
        preds.append(pred)
        X_test.append(_test_X)
        Y_test.append(_test_Y)
        
    array_pred = np.asarray(preds)
    #print '>> begin vote in overlapped patch..'
    #seg_res, possibilty_map = vote_overlapped_patch(patches_pred, index, d,h,w)
    # seconds
    elapsed = int(time.time() - start_time)

    print('Predict complete, cost [%3d] seconds' % (elapsed))

    return X_test, Y_test, array_pred


def _save_checkpoint(train_data,batch):
    td = train_data
    saver = tf.train.Saver()
    model_path = os.path.join(FLAGS.save_path,'snapshot_'+str(batch))
    
    save_path = saver.save(td.sess, model_path) #, global_step=batch)
    print("Model saved in file: %s" % save_path)   
    print ('Model saved in file: %s' % saver.last_checkpoints)

def train_model(train_data):
    td = train_data

    summaries = tf.summary.merge_all()

    # if the sess is restored from last checkpoint, do not need to 
    if FLAGS.restore_from_last:
        saver = tf.train.Saver()
        model_path = tf.train.latest_checkpoint(FLAGS.last_trained_checkpoint)
        print('training: restore last checkpoint from:%s' % model_path)
        saver.restore(td.sess, model_path)
    else:
        init_op = tf.global_variables_initializer()
        print('training: global variable initialization...')
        td.sess.run(init_op)

    lrval       = FLAGS.lr
    start_time  = time.time()
    done  = False
    batch = 0
    
    training_generator, testing_generator = get_training_and_testing_generators()

    while not done:
        batch += 1
        
        train_X, train_Y = next(training_generator)

        feed_dict = {td.X_variable : train_X, 
                     td.Y_variable : train_Y, 
                     td.learning_rate : lrval}

        cur_time = time.ctime()
        
        if batch % 10 == 0:
            ops = [td.train_minimize, td.loss, summaries,] 
            [_, loss, summary_vis] = td.sess.run(ops, feed_dict=feed_dict)

            print('[%25s], iter [%4d], Lr[%1.8f] ,loss[%3.10f]'
                            % (cur_time, batch, lrval, loss) )
            
            # Update learning rate
            if batch % FLAGS.learning_rate_reduce_life == 0:
                lrval *= FLAGS.learning_rate_percentage

            td.summary_writer.add_summary(summary_vis, batch)
            
            ############## Editted Nov 04 by Siyang Jing
            ############## Try to add validation loss
            
            val_X, val_Y = next(testing_generator)
            
            ############## Nov 14 working notes
            ## Severe bug
            ## How to enable different batch size???
            ## and different seq_length???
            
            val_feed_dict = {td.X_variable : val_X, 
                             td.Y_variable : val_Y}
            
            val_ops = [td.loss, summaries,] 
            [val_loss, val_summary] = td.sess.run(val_ops, feed_dict=val_feed_dict)

            print('[%25s], validation: iter [%4d], loss[%3.10f]'
                            % (cur_time, batch, val_loss) )
            
            td.val_sum_writer.add_summary(val_summary, batch)
            
        else:
            ops = [td.train_minimize, td.loss] 
            [_, loss] = td.sess.run(ops, feed_dict=feed_dict)
            
        if batch % FLAGS.checkpoint_period == 0:
            _save_checkpoint(td, batch)

        if batch  >= FLAGS.epoch_size:
            done = True

    _save_checkpoint(td, batch)
    print('Finished training!')
    
def prepare_dirs(delete_train_dir=False):
    ## Theoretically, should do nothing. It's just test.
    pass
    '''    
    # Create checkpoint dir (do not delete anything)
    if not tf.gfile.Exists(FLAGS.save_path):
        tf.gfile.MakeDirs(FLAGS.save_path)
    
    # Cleanup train dir
    if delete_train_dir:
        if tf.gfile.Exists(FLAGS.save_path):
            tf.gfile.DeleteRecursively(FLAGS.save_path)
        tf.gfile.MakeDirs(FLAGS.save_path)
    '''

def setup_tensorflow():
    
    config = tf.ConfigProto()
    sess = tf.Session(config=config)

    # Initialize rng with a deterministic seed
    with sess.graph.as_default():
        tf.set_random_seed(FLAGS.random_seed)
        
    random.seed(FLAGS.random_seed)
    np.random.seed(FLAGS.random_seed)

    tf.gfile.MkDir('%s/test_log' % (FLAGS.save_path,))
    summary_writer = tf.summary.FileWriter('%s/test_log' % (FLAGS.save_path,), sess.graph)
    return sess, summary_writer

class TestData(object):
    def __init__(self, dictionary):
        self.__dict__.update(dictionary)

def test():
    prepare_dirs(delete_train_dir=False)
    sess, summary_writer = setup_tensorflow()

    (X_variable, Y_variable,
     pred, loss, final_loss,
     gene_vars) = create_model(train_phase=True)

    saver = tf.train.Saver()
    ## Load from the model save path instead of the last_trained_checkpoint, 
    ## which is supposed to be the save_path of last training
    model_path = tf.train.latest_checkpoint(FLAGS.save_path)
    print('saver restore from:%s' % model_path)
    saver.restore(sess, model_path)
    
    test_data = TestData(locals())
    X_test, Y_test, array_pred = predict(test_data)
    return X_test, Y_test, array_pred


def prepare_dirs(delete_train_dir=False):
    # Create checkpoint dir (do not delete anything)
    if not tf.gfile.Exists(FLAGS.save_path):
        tf.gfile.MakeDirs(FLAGS.save_path)
    
    # Cleanup train dir
    if delete_train_dir:
        if tf.gfile.Exists(FLAGS.save_path):
            tf.gfile.DeleteRecursively(FLAGS.save_path)
        tf.gfile.MakeDirs(FLAGS.save_path)

def setup_tensorflow():
    
    config = tf.ConfigProto()
    sess = tf.Session(config=config)

    # Initialize rng with a deterministic seed
    with sess.graph.as_default():
        tf.set_random_seed(FLAGS.random_seed)
        
    random.seed(FLAGS.random_seed)
    np.random.seed(FLAGS.random_seed)

    ## Editted by Siyang Jing on Nov 4
    ## Try to add validation summary writer
    tf.gfile.MkDir('%s/training_log' % (FLAGS.save_path,))
    tf.gfile.MkDir('%s/validation_log' % (FLAGS.save_path,))
    summary_writer = tf.summary.FileWriter('%s/training_log' % (FLAGS.save_path,), sess.graph)
    val_sum_writer = tf.summary.FileWriter('%s/validation_log' % (FLAGS.save_path,), sess.graph)

    return sess, summary_writer, val_sum_writer

class TrainData(object):
    def __init__(self, dictionary):
        self.__dict__.update(dictionary)


def train():
    prepare_dirs(delete_train_dir=False)
    sess, summary_writer, val_sum_writer = setup_tensorflow()


    (X_variable, Y_variable,
     pred, loss, final_loss,
     gene_vars) = create_model(train_phase=True)

    train_minimize, learning_rate, global_step = create_optimizers(final_loss)

    train_data = TrainData(locals())
    train_model(train_data)


def main(argv=None):
    train()

if __name__ == '__main__':
    tf.app.run()