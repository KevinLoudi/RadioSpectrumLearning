# -*- coding: utf-8 -*-
"""
Created on Sun Dec 04 16:21:05 2016
Propose: Perfom build-in methods for real-data
@author: kevin
Environment: Python 2.7
"""
#############################Head##############################################
import os
dir = os.path.dirname(__file__)
#load processing modules, don't forget to rerun the load-script after modules updated
pathLoader = os.path.join(dir, 'DataLoader.py')
print os.path.exists(pathLoader)
import DataLoader as dl
pathLoader = os.path.join(dir, 'DataProcess.py')
import DataProcess as dp

#instance of module build-in class
loader=dl.DataLoader()
[time,level]=loader.loadmat()

#load data from mat files
process=dp.DataProcess(start_freq,stop_freq,step_freq,total_time)
decidecs=dp.ChannelStatus(1.0,0.0,0.5)

#############################Body##############################################

#gobal variables
start_freq=1710
stop_freq=1740
step_freq=0.025
total_time=1961

#cut the data
level_part=process.cut_by_freq(level,start_freq,stop_freq)
level_part=process.cut_by_time(level_part,0,200)

#cut the timestamps
time_part=time[0:200]

#normalize data to range [0,1]
level_part=process.normaliz_in_range(level_part,0,1)

#parse time string
time_num=process.parse_date_string(time_part)

#plot data
import numpy as np
loader.plot_spectrum_by_time(np.arange(len(time_num)),level_part[:,0])
freq=np.arange(start_freq,stop_freq,step_freq)
loader.plot_spectrum_by_freq(np.arange(len(freq)),level_part[0,:])
loader.plot_waterfall(level_part) 

##calculate a threshold decideing channel status
level_discussed=level_part[0:50,:]
loader.plot_hist(level_discussed)
threshold=decidecs.ostu_set_threshold(level_discussed)


