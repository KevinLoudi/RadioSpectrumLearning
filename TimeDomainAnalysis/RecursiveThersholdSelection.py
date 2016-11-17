# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:20:35 2016
Propose: Recursive Thresholding techiques for spectrum signal/noise decision
@author: kevin
Environment: Python 2.7
"""

import matplotlib.pyplot as plt
from skimage import data
try:
    from skimage import filters
except ImportError:
    from skimage import filter as filters
from skimage import exposure

#load data
#load level data
#current data 1710-1730MHz, step by 25kHz
import scipy.io as sio
import numpy as np
t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"].astype('float')
#camera = data.camera()
#val = filters.threshold_otsu(camera)

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

starttime=0
stoptime=50
level_part=select_part_by_freq(level,startf,stepf,stopf,1710,1730,"col")

#a sliding window that will walk through all elements of level,
#and each return the window-inside part
def sliding_window_otsu(level,win_rowsize,win_colsize,move_step):
    [rownum,colnum]=level.shape
    channelstatus=np.ndarray((rownum,colnum),float)*0
    cs_inside=np.ndarray((win_rowsize,win_colsize),float)*0
    if(win_rowsize>=rownum)and(win_colsize>=colnum)and(move_step>=min(rownum,colnum)):
        print("too large step")
        return -1
    win_ix=0
    #while(str_rowix<rownum)and(str_colix<colnum):
    for str_rowix in np.arange(0,rownum,move_step): #start,stop,step
        for str_colix in np.arange(0,colnum,move_step):
            win_inside=level[str_rowix:(str_rowix+win_rowsize),str_colix:(str_colix+win_rowsize)]
            #calculate threshold for each window
            threshold_inside=filters.threshold_otsu(win_inside)
            cs_inside=(win_inside>threshold_inside)*1
            channelstatus[str_rowix:(str_rowix+win_rowsize),str_colix:(str_colix+win_rowsize)]+=cs_inside;
            win_ix+=1
    #make a decision, under how many windows' support can an element regard 
    #as signal
    channelstatus=(channelstatus>1)*1
    print("create",win_ix," temporal windows") 
    return channelstatus              
    
#calculation
#cs->val
cs = sliding_window_otsu(level_part,20,20,5)
plt.imshow(cs)
plt,savefig('cs')
#fig=plt.imshow(cs,cmap='gray',interpolation='nearest')
#fig.show()
#fig.savefig('cs.png')

#val = filters.threshold_otsu(cs)
#hist, bins_center = exposure.histogram(cs)


#plot 
#plt.figure(figsize=(9, 4))
#plt.subplot(131)
#plt.imshow(cs, cmap='gray', interpolation='nearest')
#plt.axis('off')
#plt.subplot(132)
#plt.imshow(cs < val, cmap='gray', interpolation='nearest')
#plt.axis('off')
#plt.subplot(133)
#plt.plot(bins_center, hist, lw=2)
#plt.axvline(val, color='k', ls='--')
#
#plt.tight_layout()
#plt.show()