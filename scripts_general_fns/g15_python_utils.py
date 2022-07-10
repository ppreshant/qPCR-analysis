# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 16:41:15 2022

@author: Prashant
Purpose: Packages for python loading, commandline usage
"""

# pipes using sspipes package : simple smart pipes
# https://pypi.org/project/sspipe/
# usage: x | p(fn) | p(fn2, arg1, arg2 = px)
#--------------------------------
from sspipe import p, px # pipes

import numpy as np # numpy
import os # for directory navigation

#-----------------------------------
# change to previous directory : for going to /qPCR 
os.chdir('..')
curdir = os.getcwd()
# ------------------------------------
# dummy data for troubleshootings

# Dummy 2d array
a = np.array([[1, 2, 3], [4, 5, 6]], np.int32)

