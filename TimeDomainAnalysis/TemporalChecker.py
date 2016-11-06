# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 09:45:18 2016
Propse: tmporal folder for code validation
@author: kevin
Environment: Python 2.7
"""

#describe the inner workings of backpropagation
import numpy as np
#A neural network trained with backpropagation is attempting to use input to predict output.
#Inputs   	Output 
#0	0	1	0
#1	1	1	1
#1	0	1	1
#0	1	1	0
#We could solve this problem by simply measuring statistics between the input values and the output values.
#
#Simplest way
#X = np.array([ [0,0,1],[0,1,1],[1,0,1],[1,1,1] ])
#y = np.array([[0,1,1,0]]).T
#syn0 = 2*np.random.random((3,4)) - 1
#syn1 = 2*np.random.random((4,1)) - 1
#for j in xrange(60000):
#    l1 = 1/(1+np.exp(-(np.dot(X,syn0))))
#    l2 = 1/(1+np.exp(-(np.dot(l1,syn1))))
#    l2_delta = (y - l2)*(l2*(1-l2))
#    l1_delta = l2_delta.dot(syn1.T) * (l1 * (1-l1))
#    syn1 += l1.T.dot(l2_delta)
#    syn0 += X.T.dot(l1_delta)

# Backpropagation, in its simplest form, measures statistics like this to make a model.
# sigmoid function
def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))
    
# input dataset
X = np.array([  [0,0,1],
                [0,1,1],
                [1,0,1],
                [1,1,1] ])
    
# output dataset            
y = np.array([[0,0,1,1]]).T

# seed random numbers to make calculation
# deterministic (just a good practice)
np.random.seed(1)

# initialize weights randomly with mean 0
# first layer of weigths, connecting l0 and l1
syn0 = 2*np.random.random((3,1)) - 1

for iter in xrange(10000):

    # forward propagation
    l0 = X #first layer of the network
    l1 = nonlin(np.dot(l0,syn0)) #martix-matrix multiplication

    # how much did we miss?
    l1_error = y - l1 #second layer of the network

    # multiply how much we missed by the 
    # slope of the sigmoid at the values in l1
    l1_delta = l1_error * nonlin(l1,True) #most of the secret lay in here

    # update weights
    syn0 += np.dot(l0.T,l1_delta)

print "Output After Training:"
print l1


    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Fit Data into Distributions
#import matplotlib.pyplot as plt
#import scipy
#import scipy.stats
#size = 30000
#x = scipy.arange(size)
#y = scipy.int_(scipy.round_(scipy.stats.vonmises.rvs(5,size=size)*47))
#h = plt.hist(y, bins=range(48), color='w')
#
#dist_names = ['gamma', 'beta', 'rayleigh', 'norm', 'pareto']
#
#for dist_name in dist_names:
#    dist = getattr(scipy.stats, dist_name)
#    param = dist.fit(y)
#    pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1]) * size
#    plt.plot(pdf_fitted, label=dist_name)
#    plt.xlim(0,47)
#    print("Fit finished the ",dist_name)
#plt.legend(loc='upper right')
#plt.show()