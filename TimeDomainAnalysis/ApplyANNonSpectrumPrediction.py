# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 10:52:27 2016
Propose: Apply ANN to spectrum hole prediction
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
    occ=np.zeros([times,1])
    for i in np.arange(times):
      occ[i,0]=sum(cs[i,:])/freqs
    return occ
    
occ=calc_occ(cs)

#one-step predictor based on BP-ANN
#alldata is expected to be a one-line array,
#steps indicate how many steps is used in trainning
def rolling_window(seq, window_size):
    it = iter(seq)
    win = [it.next() for cnt in xrange(window_size)] # First window
    yield win
    for e in it: # Subsequent windows
        win[:-1] = win[1:]
        win[-1] = e
        yield win

in_set=[]
occ_seq=np.concatenate(occ)
for w in rolling_window(occ_seq, 5):
    in_set.append(w)
in_set.pop() #del in_set[-1] #[1956,1]
in_set=np.asarray(in_set)
out_set=occ[5:,].T

def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))
    
def BP_ANN_predictor(x,y):
    np.random.seed(1)
    syn0 = 2*np.random.random((5,1)) - 1

    for iter in xrange(10000):
        l0 = x #first layer of the network
        l1 = nonlin(np.dot(l0,syn0)) #martix-matrix multiplication
        l1_error = y - l1 #second layer of the network
        l1_delta = l1_error * nonlin(l1,True) #most of the secret lay in here
        syn0 += np.dot(l0.T,l1_delta)
    print "Output After Training:"
    print l1
    return

BP_ANN_predictor(in_set,out_set)
    

