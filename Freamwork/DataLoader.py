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
        

loader=DataLoader()  
print loader.count 
[time,level]=loader.loadmat()   


