# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 10:52:27 2016
Propose: Apply ANN on spectrum power prediction
@author: kevin
Environment: Python 2.7
"""

#load level data
#current data 1710-1730MHz, step by 25kHz
import scipy.io as sio
import numpy as np
t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"].astype('float')

#select part of the data by freq, presume da as a 2-d matrix
#col,row
def select_part_by_freq(da,start,step,stop,low,high,col_or_row):
    if((low<start)or(high>stop)):
        return
    lowix=int((low-start)/step)
    highix=int((high-start)/step)
    if col_or_row=="col":
        return da[:,lowix:highix]
    else:
        return da[lowix:highix,:]

#get part of level data
startf=1710
stopf=1740
stepf=0.025
level_part=select_part_by_freq(level,startf,stepf,stopf,1710,1730,"col")

#scale the loaded data into range [0,1]
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
#give an explicit range
min_max_scaler.feature_range=(0,1)
level_scaled_part = min_max_scaler.fit_transform(level_part)

#calculate channel state and show
#import matplotlib.pyplot as plt
#import pylab as pl
thre=0.35
cs=(level_scaled_part>thre)
#im = plt.matshow(cs, cmap=pl.cm.hot, aspect='auto')
#plt.colorbar(im)
#plt.show()

#calculate time-occupy rate
def calc_occ(cs):
    [times,freqs]=cs.shape
    occ=np.arange(times)*0
    for i in np.arange(times):
      occ[i]=sum(cs[i,:])/freqs
    return occ
    
occ=calc_occ(cs)

#one-step predictor based on BP-ANN
#alldata is expected to be a one-line array,
#steps indicate how many steps is used in trainning
def BP_ANN_predictor(alldata,trainpart_rate=1,steps=5):
    pre=0
    leng=alldata.size
    #shape data
    train_leng=int(leng*trainpart_rate)
#    train_set=np.zeros(((train_leng-steps),steps))
#    test_set=np.zeros((train_leng,1))
    train_set=[]
    test_set=[]

    train_ix_str=0
    train_ix_end=steps
    row_ix=0
    while(train_ix_end<train_leng-1):
#        train_set[:,row_ix]=alldata[train_ix_str:train_ix_end,1].transpose()
#        test_set[:,row_ix]=alldata[train_ix_end+1,1].transpose()
        train_set.append(alldata[train_ix_str:train_ix_end-1,0].transpose())
        test_set.append(alldata[train_ix_end,0].transpose())
        #shift 1 by row 
        train_ix_str=train_ix_str+1
        train_ix_end=train_ix_end+1
        row_ix=row_ix+1
    return pre#[train_set,test_set]#pre
    
a=BP_ANN_predictor(occ)
