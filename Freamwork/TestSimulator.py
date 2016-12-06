# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 15:48:46 2016
Propose: Perfom build-in methods for simulated-data
@author: kevin
Environment: Python 2.7
"""
import os
dir = os.path.dirname(__file__)
pathLoader = os.path.join(dir, 'Simulator.py')
import Simulator as sim

exam=sim.Simulator(40,4,10,2,0.5,10)
exam.Generate_random_signal_with_random_hole(20,100)

import scipy.special as sp
def qfunc(x):
  return 0.5 * sp.erfc(x/np.sqrt(2.0))

def qfuncinv(x):
  return np.sqrt(2) * sp.erfcinv(2*x)
  
qfunc(0.1)

