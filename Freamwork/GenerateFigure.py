# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 09:52:18 2016
Propose: Conduct ploting work
@author: kevin
Environment: Python 2.7
"""
class GenerateFigure(object):
    
    fignum=1
    fig=0
    
    def __init__(self,fignum=1):
        self.fignum=fignum
        
    def print_info(self):
        print(fignum,"res")
        return
        
        
import numpy as np
xx = np.array([-0.51, 51.2])
yy = np.array([0.33, 51.6])
means = [xx.mean(), yy.mean()]  
stds = [xx.std() / 3, yy.std() / 3]
corr = 0.8         # correlation
covs = [[stds[0]**2          , stds[0]*stds[1]*corr], 
        [stds[0]*stds[1]*corr,           stds[1]**2]] 

m = np.random.multivariate_normal(means, covs, 1000).T

x=m[0,:]
y=m[1,:]
plt.plot(x)
plt.plot(y)
np.corrcoef(x,y)

import matplotlib.pyplot as plt
plt.scatter(m[0], m[1])

        
        
