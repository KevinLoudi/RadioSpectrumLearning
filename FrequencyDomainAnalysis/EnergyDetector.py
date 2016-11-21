# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:59:56 2016
Propose: Energy detection implementation
@author: kevin
Environment: Python 2.7
"""
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

#Q-function and inverser Q function
#qfuncinv = lambda x: 0.5-0.5*sp.erf(x/np.sqrt(2))

def qfunc(x):
   return 0.5 * sp.erfc(x/np.sqrt(2.0))

def qfuncinv(x):
   return np.sqrt(2) * sp.erfcinv(2*x)

#setups
L=1000 #observe window length
snr_dB=-10 
snr=pow(10,(snr_dB/10)) #signal strength in dB and linear form
Pf=np.linspace(0.01,1,100) #pre-defined false-alarm probability
Pd=Pf*0 #detection probability
test=qfuncinv(Pf) #valid Q-function
thresh=np.ndarray((L,1),float)*0 #empty thersh group

#simulation of decter probability
for m in np.arange(Pf.size):
    print ("False Alarm Probablity set as: ")
    print (Pf[m])
    i=0
    #
    for ix in np.arange(10000):
        #noise(normal dis)+signal(normal but strong)
        n=np.random.rand(L)
        s=np.dot(np.sqrt(snr),np.random.rand(L))
        y=s+n
        
        #observed energy power
        energy=np.abs(y)**2   #pow(np.abs(y),2)
        #window-averaged energy power
        energy_ave=(1/L)*np.sum(energy)
        #np.dot((1/L),np.sum(energy))
        
        #threshold under Pf/Pd constraint
        thresh[m]=qfuncinv(Pf[m])/np.sqrt(L)+1
    
        #make decision on noise/signal
        if(energy_ave>=thresh[m]):
            i+=1  
    Pd[m]=i/float(ix);

plt.plot(Pf,Pd)

thresh=(qfuncinv(Pf)/np.sqrt(L))+1
Pd_the=qfunc((np.dot((thresh-(snr+1)),np.sqrt(L)))/(np.dot(np.sqrt(2),(snr+1))))
plt.plot(Pf,Pd_the)
        