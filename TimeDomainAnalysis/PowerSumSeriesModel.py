# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 15:54:57 2016
Propose: Updated Spectrum Timeseries Analysis, mainly on power-sum
Data:1710-1740MHz
@author: kevin
"""
import scipy.io as sio

t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"].astype('float')

#fast channel-state decided via threshold
import numpy as np
def calc_cs_via_threshold(level):
    [times,freqs]=level.shape
    cs=np.array(level.shape)
    #vital to remember to shape for row-wise max-find solution
    tmpl=np.reshape(level,(times*freqs,1))
    MaxL=max(tmpl)
    MinL=min(tmpl)
    #levels that 40% larger than the min-max gap is considered signal
    thr=MinL+(MaxL-MinL)*0.35
#    for i in np.arange(times):
#        for j in np.arange(freqs):
#            if level[i,j]>thr:
#                cs[i,j]=1
    cs=(level>thr)
    return cs,thr
    
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
    
##transform dBuV to dBm as signal strength -107
#def trans_dBuV_to_dBm(level):
#    #dBm=dBuV-107
#    #return (power=level-107)
#    
##Swept Spectrogram
#def plot_swept_spectrogram(power):
#    import pylab as pl
#    im = plt.matshow(power, cmap=pl.cm.hot, aspect='auto')
#    return
    
import pylab as pl
import matplotlib.pyplot as plt
startf=1710
stopf=1740
stepf=0.025
level_part=select_part_by_freq(level,startf,stepf,stopf,1710,1730,"col")

[cs,thr]=calc_cs_via_threshold(level_part)
im = plt.matshow(cs, cmap=pl.cm.hot, aspect='auto')
plt.colorbar(im)
plt.show()

#calculate time-occupy rate
def calc_occ(cs):
    [times,freqs]=cs.shape
    occ=np.arange(times)*0
    for i in np.arange(times):
      occ[i]=sum(cs[i,:])
    return occ
    
#sum the power of a time slot
def sum_level_by_time(level):
    [times,freqs]=level.shape
    level_sum=[]
    for i in np.arange(times):
        level_sum.append(sum(level[i,:]/100))
    return level_sum
    
occ=calc_occ(cs)
from sklearn import preprocessing
occ_scaled=preprocessing.scale(occ)
plt.plot(occ_scaled)

#ARMA model
import statsmodels.api as sm
theData=np.array(sum_level_by_time(level_part))
plt.plot(theData)
fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(311)
fig = sm.graphics.tsa.plot_acf(theData.squeeze(), lags=400, ax=ax1)
ax2 = fig.add_subplot(312)
fig = sm.graphics.tsa.plot_pacf(theData, lags=400, ax=ax2)

arma_mod30 = sm.tsa.ARMA(theData, (4,2)).fit()
print(arma_mod30.aic, arma_mod30.bic, arma_mod30.hqic)
#test if the model obey the rules
sm.stats.durbin_watson(arma_mod30.resid) #reference;;;1.9564810325209219

#observe the godness of the fit
from scipy import stats
from statsmodels.graphics.api import qqplot
resid = arma_mod30.resid
plt.plot(resid)
#test if the resid obey normal distribution
print stats.normaltest(resid) #reference (49.845005822356519, 1.5007021697683078e-11)
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
fig = qqplot(resid, line='q', ax=ax, fit=True)

#to see if resid have any time-correlation
fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
fig = sm.graphics.tsa.plot_acf(resid.squeeze(), lags=40, ax=ax1)
ax2 = fig.add_subplot(212)
fig = sm.graphics.tsa.plot_pacf(resid, lags=40, ax=ax2)

predict_sunspots = arma_mod30.predict(dynamic=True)
print(predict_sunspots)
