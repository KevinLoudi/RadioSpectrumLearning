# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 15:12:04 2016
Conduct simulation on signal-generating, detecting
@author: kevin
Environment: Python 3.5
"""

class Simulator(object):
    
    signal_mean=0.0
    signal_std=0.0
    noise_mean=0.0
    noise_std=0.0
    vancy_rate=0.5
    signal_noise_rate=10 #dB
    
    def __init__(self,signal_mean=0.0,signal_std=0.0,noise_mean=0.0,noise_std=0.0,vancy_rate=0.5,signal_noise_rate=10):
        self.signal_mean=signal_mean
        self.signal_std=signal_std
        self.noise_mean=noise_mean
        self.noise_std=noise_std
        self.vancy_rate=vancy_rate
        self.signal_noise_rate=signal_noise_rate;
             
    #signal simulation continuing 'times' time, each time with 'slot' slots 
    def Generate_random_signal_with_random_hole(self,times,slot):
        #as a build modules
        n_mean=self.noise_mean
        n_std=self.noise_std
        s_mean=self.signal_mean
        s_std=self.signal_std
        vac_rate=self.vancy_rate
        snr=self.signal_noise_rate
        
        #as a test block
        times=20
        slot=100
        n_mean=10.0 
        n_std=2.0
        s_mean=40.0
        s_std=5.0
        vac_rate=0.5
        snr=-10.0
        
        #L: total time length // Num: visit times
        #model the random access from primary user
        def random_visit(L,Num):
            #random inital signal start time
            start_time=np.random.uniform(0,L,Num)
            start_time=np.sort(start_time,kind='mergesort',axis=0)
            #intend to preserve the longest possible stay time of the signal's visit
            max_possible_gap=np.zeros(len(start_time))
            #signal visit time
            stay_time=np.zeros(len(max_possible_gap))
            endix=len(start_time)-1
            
            for i in np.arange(0,len(start_time)-1,1):
                max_possible_gap[i]=int(start_time[i+1]-start_time[i])
                #random set visit time within range [1,max_possibe]
                stay_time[i]=int(np.random.uniform(1,max_possible_gap[i],1))
            #special handle for the end item
            max_possible_gap[endix]=int(L-start_time[endix])
            stay_time[endix]=int(np.random.uniform(1,max_possible_gap[endix],1))
        
            #get level time of each visit
            level_time=start_time+stay_time
            rcs=np.zeros(L)
            
            #set the status of signal-presence as '1'
            for ix in np.arange(len(start_time)):
                for i in np.arange(L):
                    if i>start_time[ix] and i<level_time[ix]:
                        rcs[i]+=1

            plt.plot(rcs)
            return rcs
            
        random_visit(10000,3)
        
        def Multi_signal_generator():
            power_arr=np.arange(15,45,3)    
            s_list=[]
            rcs_list=[]
            L=10000
            #signal source 1
            ix=0
            rcs_list.append(random_visit(L,10))
            s=rcs_list[ix]*(np.random.normal(power_arr[ix],power_arr[ix]/10,L))
            
            
            
            
            start_time=start_time.tolist()
            start_time.
            return 0
        
        import numpy as np
        import matplotlib.pyplot as plt
        
        #snr_lin=pow(10,(snr/10)) #signal strength in dB and linear form
        Y=np.zeros(slot*times) #Y=S+N
        S=np.zeros(slot*times)
        N=np.zeros(slot*times)
        RCS=np.zeros(slot*times) #realy channel state
        
        #M-K simulation
        for t in np.arange(times):
            cs=(np.random.rand()>vac_rate) #traffic condition
            if(cs):
                #H1
                #s=np.dot(np.sqrt(snr_lin),np.random.normal(s_mean,s_std,slot))
                s=np.random.normal(s_mean,s_std,slot)
                n=np.random.normal(n_mean,n_std,slot)
                Y[t*slot:(t+1)*slot]=s+n
                S[t*slot:(t+1)*slot]=s
                N[t*slot:(t+1)*slot]=n
                RCS[t*slot:(t+1)*slot]=np.ones(slot)   
            else:
                #H1
                n=np.random.normal(n_mean,n_std,slot)
                Y[t*slot:(t+1)*slot]=n
                S[t*slot:(t+1)*slot]=np.zeros(slot)
                N[t*slot:(t+1)*slot]=n
                RCS[t*slot:(t+1)*slot]=np.zeros(slot) 
            
        plt.figure(1)
        plt.plot(N,'g')
        plt.plot(Y,'r')
        #plt.plot(RCS,'k')
        plt.show()
            
        plt.figure(2)
        plt.hist(Y)
        
        try:
            from skimage import filters
        except ImportError:
            from skimage import filter as filters
        thre=filters.threshold_otsu(Y)
        print("threshold_ostu: ",thre)
        
        #false alarm probability
        import scipy.special as sp
        def qfunc(x):
            return 0.5 * sp.erfc(x/np.sqrt(2.0))

        def qfuncinv(x):
            return np.sqrt(2) * sp.erfcinv(2*x)
        
        #variation of detection performance with changes on observe window
        Num_vec=np.arange(0.01,slot,0.02,dtype='float')
        Pr_detect=np.arange(0.01,slot,0.02,dtype='float')*0
        Pr_false_alarm=np.arange(0.01,slot,0.02,dtype='float')*0
        ix=0
        for Num in Num_vec:
            tmp=((thre-Num*(s_std**2+n_std**2)))/(np.sqrt(2*Num*(s_std**2+n_std**2)**2))
            Pr_detect[ix]=qfunc(tmp)
            Pr_false_alarm[ix]=qfunc((thre-Num*n_std**2)/(np.sqrt(2*Num*n_std**4)))
            ix+=1
        fig = plt.figure()
        plt.plot(Pr_false_alarm,Pr_detect)
        ax = fig.add_subplot(111)
        ax.axis([0,1,0,1])
            
#        Num=1 #sensing period
#        tmp=((thre-Num*(s_std**2+n_std**2)))/(np.sqrt(2*Num*(s_std**2+n_std**2)**2))
#        print(tmp)
#        Pr_detect=qfunc(tmp)
#        Pr_false_alarm=qfunc((thre-Num*n_std**2)/(np.sqrt(2*Num*n_std**4)))
#        print("detection probability: ",Pr_detect)
#        print("false alarm probability: ",Pr_false_alarm)
        
        #classify Y into signal and noise
        #estimate features of signal and noise
        Y_lst=Y.tolist()
        Y_only_noise=filter(lambda x: x<thre, Y_lst)
        Y_with_signal=filter(lambda x: x>=thre, Y_lst)
        Y_noise_mean=np.mean(Y_only_noise)
        Y_noise_std=np.sqrt(np.var(Y_only_noise))
        Y_with_signal=np.asarray(Y_with_signal)
        Y_with_signal-=np.ones(len(Y_with_signal))*Y_noise_mean
        Y_signal_mean=np.mean(Y_with_signal)
        Y_signal_std=np.sqrt(np.var(Y_with_signal))
        
        
        import scipy.signal as ss
        import matplotlib.pyplot as plt
        t = np.linspace(0, 1, 500, endpoint=False)
        plt.plot(t, signal.square(2 * np.pi * 5 * t))
        plt.ylim(-2, 2)
        
        return 0
            
            
        

