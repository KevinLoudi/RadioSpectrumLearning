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
            
        #Num:source number, L: total time length
        def Multi_signal_generator(Num,L,occur_time):
            power_arr=np.arange(15,45,30.0/Num)    
            s_list=[]
            rcs_list=[]
            #signal source 1
            ix=0
            while(ix<Num):
                rcs_list.append(random_visit(L,occur_time))
                s_list.append(rcs_list[ix]*(np.random.normal(power_arr[ix],power_arr[ix]/10,L)))
                #plt.plot(s)
                ix+=1
            return s_list
            
            
        #generate uniform noise
        def uniform_noise_generator(noise_mean,noise_std,L):
            import numpy as np
            n=np.random.normal(noise_mean,noise_std,L)
            return n    
            
        #use recursive threshold
        def recursive_oneside_hypthesis_testing(data_vec,max_iter,z_alph):    
               #z_alph=0.6#1.645 #95% confidence
               d_mean=list()
               d_std=list()
               d_cutpoint=list()
               for ix in np.arange(max_iter):
                   #collect mean and std of each iterate
                   d_mean.append(np.mean(data_vec))
                   d_std.append(np.std(data_vec))
                   d_cutpoint.append(d_mean[ix]+z_alph*d_std[ix])
                   #subtract the signal part with 95% confidence
                   #and re-reshape to array
                   data_list=data_vec.tolist()
                   data_list_filter=filter(lambda x: x<d_cutpoint[ix], data_list)
                   data_vec=np.asarray(data_list_filter)
               print ix
               print d_cutpoint[ix]
               return [d_std,d_mean,d_cutpoint,ix]

        #tail function for normal distribution
        def qfunc(x):
              import scipy.special as sp
              return 0.5 * sp.erfc(x/np.sqrt(2.0))

        def qfuncinv(x):
             import scipy.special as sp
             return np.sqrt(2) * sp.erfcinv(2*x)
             
        #calculate detect probability needs threshold, signal std
        #noise std and observe window size
        def calc_p_d(thre,s_std,n_std,obs_size):
            import numpy as np
            return qfunc((thre-obs_size*(s_std**2+n_std**2))/(np.sqrt(2*obs_size*(s_std**2+n_std**2)**2)))
        
        #calculate false alarm prrobability needs threshold, noise std 
        #and observe window size
        def calc_p_fa(thre,n_std,obs_size):
            import numpy as np
            return qfunc((thre-obs_size*n_std**2)/(np.sqrt(2*obs_size*n_std**4)))
        
        #generate a random array y, y is related to x and have a correlation-rate of r
        def correlated_value(x,r):
            r2=r**2
            var=1-r2
            std=np.sqrt(var)
            err=np.random.normal(0.0,std,len(x))
            y=r*x+err
            return y
            
        #freqency length, time span, visit time of each signal source, time corre,noise parameters
        def generate_spectrum(freqlen=1000,timelen=10,visit_time=3,time_corr=0.98,noise_mean=10,noise_std=2):
            L=freqlen #length of each fream of spectrum
            T=timelen #signal source number
            s=Multi_signal_generator(T,L,visit_time)
            ss=np.zeros(L)
            for ix in np.arange(len(s)):
                ss=ss+s[ix]
            n=uniform_noise_generator(noise_mean,noise_std,L)
            y=ss+n #signal and noise
            yt=np.zeros([T,L])
            ix=0
            tmp=ss
            while(ix<T):
                yt[ix,:]=tmp
                tmp=correlated_value(tmp,time_corr)
                ix+=1
            import matplotlib.pyplot as plt
            import pylab as pl
            im = plt.matshow(yt, cmap=pl.cm.hot, aspect='auto')
            plt.colorbar(im)
            plt.show()
            return yt
             
    #signal simulation continuing 'times' time, each time with 'slot' slots 
       def Generate_random_signal_with_random_hole(self,times,slot):
        #as a build modules
        n_mean=self.noise_mean
        n_std=self.noise_std
        s_mean=self.signal_mean
        s_std=self.signal_std
        vac_rate=self.vancy_rate
        snr=self.signal_noise_rate
        
        #simulate with 3 sourc in a length of 10000
        #as a test block
        import numpy as np
        import matplotlib.pyplot as plt
        #parameters for signal and noise
        times=20
        slot=100
        n_mean=10.0 
        n_std=2.0
        s_mean=40.0
        s_std=5.0
        vac_rate=0.5
        snr=-10.0 
        
        #generate signal with different strength
        L=10000
        Source_num=10
        random_visit(L,Source_num)
        s_list=Multi_signal_generator(Source_num,L)
        
        #combine signal from sources
        ss=np.zeros(L)
        for ix in np.arange(len(s_list)):
            ss=ss+s_list[ix]
        
        nn=uniform_noise_generator(n_mean,n_std,L)
        yy=ss+nn
        #received signal from several sources and with noise
        plt.plot(yy)
        plt.hist(yy,bins=int(max(yy)))
        y_std=np.std(yy)
        y_mean=np.mean(yy)
        #provide a decision threshold through histgram
        #use global threshold
        try:
            from skimage import filters
        except ImportError:
            from skimage import filter as filters
        thre=filters.threshold_otsu(yy)
        print("threshold_ostu: ",thre)
        cs_est=np.zeros(len(yy))
        cs_est=yy>thre
        
        t=np.arange(0,L)
        plt.plot(t,yy,'g',t,cs_est*70,'r')
        plt.plot(cs_est)

        [y_std,y_mean,y_cutpoint,ix]=recursive_oneside_hypthesis_testing(yy,20,0.8) 
        plt.plot(y_cutpoint) #cutpoint will converge to the threshold
        
        #estimate signal and noise character
        thre=y_cutpoint[19]
        obs_size=10
        ss_mean=np.mean(ss)
        ss_std=np.std(ss)
        nn_mean=np.mean(ss)
        nn_std=np.std(nn)
        p_d=calc_p_d(thre,ss_std,nn_std,obs_size)
        p_fa=calc_p_fa(thre,nn_std,obs_size)
        
           
        #relationship between P_d and P_fa when both signal and 
        #noise follows normal distribution
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
        
        #decided threshold through P_d , P_fa point view
        P_fa=0.01
        n_std_est=n_std*0.96
        obs_len=20
        thre=(qfuncinv(P_fa)*np.sqrt(2*obs_len)+obs_len)*n_std_est
        print(thre)
        
        
        ##generate another array that is related to ss        
        time_num=30
        yyy=np.zeros([time_num,len(yy)])
        ix=0
        tmp=ss
        while(ix<channel_num):
            yyy[ix,:]=tmp
            tmp=correlated_value(tmp,0.95)
            ix+=1
            
        import matplotlib.pyplot as plt
        import pylab as pl
        im = plt.matshow(yyy, cmap=pl.cm.hot, aspect='auto')
        plt.colorbar(im)
        plt.show()
        
            
#        import numpy as np
#        import matplotlib.pyplot as plt
#        
#        #snr_lin=pow(10,(snr/10)) #signal strength in dB and linear form
#        Y=np.zeros(slot*times) #Y=S+N
#        S=np.zeros(slot*times)
#        N=np.zeros(slot*times)
#        RCS=np.zeros(slot*times) #realy channel state
#        
#        #M-K simulation
#        for t in np.arange(times):
#            cs=(np.random.rand()>vac_rate) #traffic condition
#            if(cs):
#                #H1
#                #s=np.dot(np.sqrt(snr_lin),np.random.normal(s_mean,s_std,slot))
#                s=np.random.normal(s_mean,s_std,slot)
#                n=np.random.normal(n_mean,n_std,slot)
#                Y[t*slot:(t+1)*slot]=s+n
#                S[t*slot:(t+1)*slot]=s
#                N[t*slot:(t+1)*slot]=n
#                RCS[t*slot:(t+1)*slot]=np.ones(slot)   
#            else:
#                #H1
#                n=np.random.normal(n_mean,n_std,slot)
#                Y[t*slot:(t+1)*slot]=n
#                S[t*slot:(t+1)*slot]=np.zeros(slot)
#                N[t*slot:(t+1)*slot]=n
#                RCS[t*slot:(t+1)*slot]=np.zeros(slot) 
#            
#        plt.figure(1)
#        plt.plot(N,'g')
#        plt.plot(Y,'r')
#        #plt.plot(RCS,'k')
#        plt.show()
#            
#        plt.figure(2)
#        plt.hist(Y)
#        
#        try:
#            from skimage import filters
#        except ImportError:
#            from skimage import filter as filters
#        thre=filters.threshold_otsu(Y)
#        print("threshold_ostu: ",thre)
        
        #false alarm probability

            
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
            
            
        

