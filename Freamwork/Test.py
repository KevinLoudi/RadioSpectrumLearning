# -*- coding: utf-8 -*-
"""
Created on Sun Dec 04 16:21:05 2016
Propose: Perfom build-in methods
@author: kevin
Environment: Python 2.7
"""
import os
dir = os.path.dirname(__file__)
pathLoader = os.path.join(dir, 'DataLoader.py')
print os.path.exists(pathLoader)
import DataLoader as dl
pathLoader = os.path.join(dir, 'DataProcess.py')
import DataProcess as dp

loader=dl.DataLoader()
[time,level]=loader.loadmat()

process=dp.DataProcess()
level_part=process.cut_by_freq(level,88,108)
