# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:07:06 2016

@author: kevin
"""


def works_fine():

    a = 5

    b = 6

    assert(a + b == 11)



def throws_an_exception():

    a = 5

    b = 6

    assert(a + b == 10)



def calling_things():

    works_fine()

    throws_an_exception()



calling_things()
