import os
import sys
import tensorflow as tf
import configparser

import matlab.engine

MATLAB_ENGINE = matlab.engine.start_matlab()
FLAGS = tf.app.flags.FLAGS

## Configuration File Parse
CONFIG_DIR = './config.ini'
if len(sys.argv)>2 and sys.argv[2][-4:]=='.ini':
    CONFIG_DIR = sys.argv[2]
CFP = configparser.ConfigParser()
CFP.read(CONFIG_DIR)


def ResetValue(section,name,value):
    CFP[section].update({name:value})
    
    
def ReadConfigFile(path=CONFIG_DIR):
    CFP.read(path)