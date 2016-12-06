# -*- coding: utf-8 -*-
"""
Created on Sun Dec 04 15:58:26 2016
Propose: File I/O methods
@author: kevin
Environment: Python 2.7
"""

class DataLoader(object):
    "This is a class for data loading method"
    
    count=0;
    
    def __init__(self,count=0):
        self.count=count
    
    ##load data from a setted folder
    def loadmat(self):
        import scipy.io as sio
        import os
        dir = os.path.dirname(__file__)
        filename_l = os.path.join(dir, 'TestData/level.mat')
        filename_t = os.path.join(dir, 'TestData/timestamp.mat')
        if(os.path.exists(filename_l)and os.path.exists(filename_t)):
            t=sio.loadmat(filename_t)
            time=t["dateStamp"]
            l=sio.loadmat(filename_l)
            level=l["dataLevel"].astype('float')
            return [time,level]
        else:
            return -1
    
    #load data in a form of figure
    def plot_spectrum_by_time(self,x,y):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        #here only plot one figure
        fig.suptitle('spectrum time', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
        ax.set_title('one')
        ax.set_xlabel('time')
        ax.set_ylabel('power')
        ax.plot(x,y)
        #set scope for demo
        ax.axis([0,len(x)+3,0,max(y)])
        plt.show()
        return 0 
        
    def plot_spectrum_by_freq(self,x,y):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        #here only plot one figure
        fig.suptitle('spectrum freq', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
        ax.set_title('two')
        ax.set_xlabel('freq')
        ax.set_ylabel('power')
        ax.plot(x,y)
        #set scope for demo
        ax.axis([0,len(x)+3,0,max(y)])
        plt.show()
        return 0 
        
    def plot_waterfall(self,data):
        import matplotlib.pyplot as plt
        import pylab as pl
        im = plt.matshow(data, cmap=pl.cm.hot, aspect='auto')
        plt.colorbar(im)
        plt.show()
        
    def plot_hist(self,data):
        import numpy as np
        import matplotlib.pyplot as plt
        data_vec=np.reshape(data,data.size)
        plt.hist(data_vec,bins='auto')
        return "histgram plotted"


