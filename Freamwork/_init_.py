# -*- coding: utf-8 -*-
"""
Created on Sun Dec 04 16:34:01 2016

@author: kevin
"""

if __name__ != '__main__':
    from .DataLoader import *  # defines __all__
    __all__.remove('__all__')  # prevent export (optional)