# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 09:15:18 2016
Propose: Fit power data into distribution
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

def select_part_by_freq(da,start,step,stop,low,high,col_or_row):
    if((low<start)or(high>stop)):
        return
    lowix=int((low-start)/step)
    highix=int((high-start)/step)
    if col_or_row=="col":
        return da[:,lowix:highix]
    else:
        return da[lowix:highix,:]

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
import matplotlib.pyplot as plt
import pylab as pl
thre=0.35
cs=(level_scaled_part>thre)
im = plt.matshow(cs, cmap=pl.cm.hot, aspect='auto')
plt.colorbar(im)
plt.show()

#calculate time-occupy rate
def calc_occ(cs):
    [times,freqs]=cs.shape
    occ=np.arange(times)*0
    for i in np.arange(times):
      occ[i]=sum(cs[i,:])/freqs
    return occ
    
#occ=calc_occ(cs)
#plt.plot(occ)

#fit power data in find the best theoretical distribution model
import matplotlib.pyplot as plt
import scipy
import scipy.stats

dist_names = ['gamma', 'beta']#,'rayleigh', 'norm', 'pareto']
#size of x and y expected to be indentical
x=np.arange(1710,1730,0.025)
y=level_part[400,:]
size=x.size()
h = plt.hist(y, bins=range(800), color='w')

for dist_name in dist_names:
    dist = getattr(scipy.stats, dist_name)
    param = dist.fit(y)
    pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1]) * size
    plt.plot(pdf_fitted, label=dist_name)
    plt.xlim(0,799)
    print("Fit finished the ",dist_name)
plt.legend(loc='upper right')
plt.show()

