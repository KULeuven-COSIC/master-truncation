#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 14:18:18 2023

@author: martinzbudila
"""
import numpy as np
from tensorflow import keras

(X_train, y_train), (X_test, y_test) = keras.datasets.mnist.load_data()

X_train = X_train.reshape(X_train.shape[0], -1)
X_test = X_test.reshape(X_test.shape[0], -1)

y_test = np.eye(y_test.shape[0],10)[y_test]
y_test = y_test.reshape(y_test.shape[0],10)

y_train = np.eye(y_train.shape[0],10)[y_train]
y_train = y_train.reshape(y_train.shape[0], 10)

X_train = X_train/255
X_test = X_test/255

def get_data_files():
    np.savetxt("Y_val.txt", np.array(y_test), comments = '')    

    np.savetxt("X_val.txt", np.array(X_test), comments = '') 

    np.savetxt("Y_train.txt", np.array(y_train), comments = '')    

    np.savetxt("X_train.txt", np.array(X_train), comments = '')  