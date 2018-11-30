import os
import sys
import tensorflow as tf
import configparser
FLAGS = tf.app.flags.FLAGS

## Configuration File Parse
CONFIG_DIR = './config.ini'
if len(sys.argv)>1 and sys.argv[1][-4:]=='.ini':
    CONFIG_DIR = sys.argv[1]
CFP = configparser.ConfigParser()
CFP.read(CONFIG_DIR)
