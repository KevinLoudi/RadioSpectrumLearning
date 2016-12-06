# -*- coding: utf-8 -*-
"""
Created on Sun Dec 04 16:21:05 2016
Propose: Perfom build-in methods
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
#import numpy as np
#loader.plot_spectrum_by_time(np.arange(len(time_num)),level_part[:,0])
#freq=np.arange(start_freq,stop_freq,step_freq)
#loader.plot_spectrum_by_freq(np.arange(len(freq)),level_part[0,:])
#loader.plot_waterfall(level_part) 

##calculate a threshold decideing channel status
#level_discussed=level_part[0:50,:]
#loader.plot_hist(level_discussed)
#threshold=decidecs.ostu_set_threshold(level_discussed)

#from energy-detection point view
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

plt.figure

Times=20
Slot=100
snr_dB=10
snr=pow(10,(snr_dB/10)) #signal strength in dB and linear form
y=np.arange(Slot*Times,dtype='float')*0
N=np.arange(Slot*Times,dtype='float')*0
RCS=np.arange(Slot*Times,dtype='float')*0
   
#test
#a=random_vec_standard(10,2,300)
#a_std=np.sqrt(np.var(a))
#a_mean=np.mean(a)

#variance of signal 4, variance of noise 2
s_std=4
s_mean=35
n_std=6
n_mean=10
for t in np.arange(Times):
    cs=(np.random.rand()>0.5) #traffic condition
    if(cs):
        #H1
        s=np.dot(np.sqrt(snr),decidecs.random_vec_standard(s_mean,s_std,Slot))
        n=decidecs.random_vec_standard(n_mean,n_std,Slot)
        y[t*Slot:(t+1)*Slot]=s+n
        N[t*Slot:(t+1)*Slot]=n
        RCS[t*Slot:(t+1)*Slot]=np.ones(Slot)
    else:
        #H0
        n=decidecs.random_vec_standard(n_mean,n_std,Slot)     
        y[t*Slot:(t+1)*Slot]=n
        N[t*Slot:(t+1)*Slot]=n
        RCS[t*Slot:(t+1)*Slot]+=np.ones(Slot)*0
    
plt.plot(N,'g')
plt.plot(y,'r')
plt.plot(RCS,'y')
plt.show()

loader.plot_hist(y)
thre=decidecs.ostu_set_threshold(y)

#false alarm probability
def qfunc(x):
   return 0.5 * sp.erfc(x/np.sqrt(2.0))

def qfuncinv(x):
   return np.sqrt(2) * sp.erfcinv(2*x)

Num=Times*Slot
Pr_detect=qfunc((thre-Num*(s_std*s_std+n_std*n_std))/(np.sqrt(2*Num*(s_std*s_std+n_std*n_std))))
Pr_false_alarm=qfunc((thre-Num*n_std*n_std)/(np.sqrt(2*Num*n_std*n_std*n_std*n_std)))

