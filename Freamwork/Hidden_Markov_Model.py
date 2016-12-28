# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 10:54:53 2016
Propose: Hidden Markov Model
@author: kevin
Environment: Python 3.5
"""
#build HMM 
import hmmlearn.hmm as hmm
import numpy as np

transmat = np.array([[0.7, 0.3],
                      [0.3, 0.7]])
emitmat = np.array([[0.9, 0.1],
                    [0.2, 0.8]])

startprob = np.array([0.5, 0.5])

h = hmm.MultinomialHMM(n_components=2, startprob=startprob,transmat=transmat)
h.emissionprob_ = emitmat
# works fine
h.fit([[0, 0, 1, 0, 0]]) 
# h.fit([[0, 0, 1, 0, 0], [0, 0], [1,1,1]]) # this is the reason for such 
                                            # syntax, you can fit to multiple
                                            # sequences    
print(h.decode([0, 0, 1, 0, 0]))
print(h)