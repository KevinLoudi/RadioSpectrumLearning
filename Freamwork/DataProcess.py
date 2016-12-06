# -*- coding: utf-8 -*-
"""
Created on Sun Dec 04 17:00:55 2016
Propose: Methods for general data process
@author: kevin
Environment: Python 2.7
"""

class DataProcess(object):
    
    start_freq=20
    stop_freq=3000
    step_freq=0.025
    time_length=100
    
    def __init__(self,start_freq=20,stop_freq=3000,step_freq=0.025,time_length=100):
        self.start_freq=start_freq
        self.stop_freq=stop_freq
        self.step_freq=step_freq
        self.time_length=time_length
        
    def cut_by_freq(self,data,low=int,high=int):
        if(low<self.start_freq or high>self.stop_freq or data.shape==(1,1)):
          print("error input:","out of bound","abnormal input data")
          return -1 #error input,out of bound and error shape
        else:
          if(low>high):
           #swap low with high
           tmp=low
           low=high
           high=tmp
          else:
           lowix=(low-self.start_freq)/self.step_freq
           highix=(high-self.start_freq)/self.step_freq
           return data[:,lowix:highix] 

    def cut_by_time(self,data,low=int,high=int):
        if(low<0 or high>self.time_length or data.shape==(1,1)):
          print("error input:","out of bound","abnormal input data")
          return -1 #error input,out of bound and error shape
        else:
          if(low>high):
           #swap low with high
           tmp=low
           low=high
           high=tmp
          else:
           return data[low:high,:]

    def normaliz_in_range(self,data,floor=int,ceil=int):
        #scale the loaded data into range [0,1]
        from sklearn import preprocessing
        min_max_scaler = preprocessing.MinMaxScaler()
        #give an explicit range
        min_max_scaler.feature_range=(floor,ceil) 
        return min_max_scaler.fit_transform(data)
    
    #parse date string to num-form date
    def parse_date_string(self,time,timestr='%d-%b-%Y %H:%M:%S'):
        from datetime import datetime
        #initialize a datetime
        time_num=[]
        ix=0
        while(ix<len(time)):
            time_num.append(datetime.strptime(time[ix],timestr)) #like '15-Dec-2015 19:00:00' 
            ix+=1
        print("Results like ",time_num[0])
        return time_num
        
        
class ChannelStatus(object):
    maxx=1.0
    minn=0.0
    threshold=0.5
    
    def __init__(self,maxx=1.0,minn=0.0,threshold=0.5):
        self.maxx=maxx
        self.minn=minn
        self.threshold=threshold
    
    def manual_set_threshold(self,threshold=0.5):
        self.threshold=threshold
        return 0

    def ostu_set_threshold(self,data):
        import numpy as np
        try:
            from skimage import filters
        except ImportError:
            from skimage import filter as filters
        
        if(not isinstance(data,np.ndarray)):
            return -999
        data_vec=np.reshape(data,data.size)
        threshold_ostu=filters.threshold_otsu(data_vec)
        print("threshold_ostu: ",threshold_ostu)
        self.threshold=threshold_ostu
        return threshold_ostu
    
    #generate random numbers in uniform distribution
    def random_vec_standard(self,mean,std,length):
        import numpy as np
        vec=np.random.normal(mean,std,length) #[low high size]
        return vec

    #generate a vector of random number in uniform distribution,
    #that is between [low,high]
    def random_vec_uniform(self,low,high,length):
        import numpy as np
        vec=np.random.uniform(low, high,length) #[low high size]
        return vec

        