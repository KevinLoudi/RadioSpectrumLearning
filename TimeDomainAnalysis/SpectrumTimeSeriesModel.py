# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 13:25:45 2016
Propose: Model and Predict Spectrum Power Time Series
with ARIMA
Enviroment: Python 2.7
@author: kevin
"""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans,vq

#provide relative path
dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'TestData/DayChannelState.csv')

#load csv file for channelstate
res=pd.read_csv(filename)
#channel state in '0' or '1'
cs=np.array(res).astype('float')

#calculate channe-state via k-means
def cluster_bytime(level):
    # computing K-Means with K = 2 (2 clusters)
    centroids,_ = kmeans(level,2)
    # assign each sample to a cluster
    idx,_ = vq(level,centroids)
    return idx #np.transpose(idx)
    
#clustering by time slice
def calc_cs(lev):
    [row,col]=lev.shape
    tlevave=np.transpose(np.arange(col)*0.0)
    tlevmax=tlevave
    tlevmin=tlevave
    for i in np.arange(row):
        #cluster act,return a classify label "channel" or "dumy"
        #do clustering for each time slot
        tmprow=cluster_bytime(lev[i,:])
        tlevave=tlevave+tmprow
        for j in np.arange(col):
            if tlevmax[j]<tmprow[j]:   
               tlevmax[j]=tmprow[j]
            if tlevmin[j]>tmprow[j]:
               tlevmin[j]=tmprow[j] 
    #Here I want to see how like the sample could represent a channel
    tlevave=tlevave/row #the possibilty of being a channel
    return tlevmin,tlevave,tlevmax

#calculate time-occupy rate
def calc_occ(cs):
    [times,freqs]=cs.shape
    occ=np.arange(times)*0
    for i in np.arange(times):
      occ[i]=sum(cs[i,:])
    return occ

#load more data
t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"]
plt.plot(level[:,100])

[lmin,lmean,lmax]=calc_cs(level)

