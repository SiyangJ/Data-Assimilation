import tensorflow as tf
from tensorflow import keras

import numpy as np

from EnKF_func import *
from AuxFuncs import *
import matplotlib.pyplot as plt

def ML(train_data,train_labels,
       test_data,test_labels,
       model=None,
       MLDEBUG=True,
       MLPLOTTING=False,
       EPOCHS=1000):
    
    # Default Model
    if model is None:
        def build_model():
            model = keras.Sequential([
            keras.layers.Dense(20, activation=tf.nn.relu, 
                               input_shape=(train_data.shape[1],)),
            keras.layers.Dense(20, activation=tf.nn.relu),
            keras.layers.Dense(20)
          ])

            #optimizer = tf.train.RMSPropOptimizer(0.001)

            model.compile(loss='mse',
                        optimizer='adam',
                        metrics=['mae','acc'])
            return model

        model = build_model()

    # Display training progress by printing a single dot for each completed epoch.
    class PrintDot(keras.callbacks.Callback):
        def on_epoch_end(self,epoch,logs):
            if epoch % 100 == 0: print('')
            print('.', end='')

    # Store training stats
    history = model.fit(train_data, train_labels, epochs=EPOCHS,
                        validation_split=0.2, verbose=0,
                        callbacks=[PrintDot()] if MLDEBUG else None)
    
    [loss, mae, acc] = model.evaluate(test_data, test_labels, verbose=0)

    if MLPLOTTING:
        def plot_history(history):
            plt.figure()
            plt.xlabel('Epoch')
            plt.ylabel('Mean Abs Error')
            plt.plot(history.epoch, np.array(history.history['mean_absolute_error']), 
                   label='Train Loss')
            plt.plot(history.epoch, np.array(history.history['val_mean_absolute_error']),
                   label = 'Val loss')
            plt.legend()
            plt.ylim([0,5])

        plot_history(history)

        def plot_acc_history(history):
            plt.figure()
            plt.xlabel('Epoch')
            plt.ylabel('Mean Rel Error')
            plt.plot(history.epoch, np.array(history.history['acc']), 
                   label='Train Loss')
            plt.plot(history.epoch, np.array(history.history['val_acc']),
                   label = 'Val loss')
            plt.legend()
            plt.ylim([0,5])

        plot_acc_history(history)

        print("")
        print("Testing set Mean Abs Error: {:4.2f}".format(mae))
        print("Testing set Accuracy: {:4.2f}".format(acc))
    
    return (model,mae,acc)

