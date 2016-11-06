# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 09:15:18 2016
Propose: Fit power data into distribution
@author: kevin
Environment: Python 2.7
"""
#load level data
import scipy.io as sio
t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"].astype('float')

#scale the loaded data into range [0,1]
from sklearn import preprocessing
import numpy as np
min_max_scaler = preprocessing.MinMaxScaler()
#give an explicit range
min_max_scaler.feature_range=(0,1)
level_scaled = min_max_scaler.fit_transform(level)