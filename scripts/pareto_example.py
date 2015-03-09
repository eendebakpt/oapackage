# -*- coding: utf-8 -*-
"""
Python performance test

Log with:
    
    > python ptest.py > test0-ddmmyyy.txt
"""

#%% Load necessary packages 
import sys
import os
import numpy as np
import time
import glob
import matplotlib.pyplot as plt
import subprocess
import time
import getopt
import platform


x=np.random.rand(2, 50)
#x=np.sqrt(x)
#x=np.square(x)

for ii in range(0, x.shape[1]):
    w=x[:,ii]
    fac=.6+.4*sqrt(w[0]**2+w[1]**2)
    x[:,ii]=(1/fac)*w

import oalib

par=oalib.ParetoDoubleLong()

for ii in range(0, x.shape[1]):
    w=oalib.doubleVector( (x[0,ii], x[1,ii]))
    par.addvalue(w, ii)

lst=par.allindices()

xp=x[:,lst]

plt.figure(1)
plt.clf()
h=plt.plot(x[0,:], x[1,:], '.b', markersize=16, label='Non Pareto-optimal')
hp=plt.plot(xp[0,:], xp[1,:], '.r', markersize=16, label='Pareto optimal')
plt.xlabel('Parameter 1', fontsize=16)
plt.ylabel('Parameter 2', fontsize=16)
plt.xticks([])
plt.yticks([])
plt.legend(loc=3, numpoints=1)
#lt.axis('off')
#plt.legend([h,hp], ['Non-pareto', 'Pareto'])

#plt.savefig('/home/eendebakpt/misc/homepage/files/pareto-example.png', transparent=True)

if 0:
    plt.figure(2)
    h=plt.plot(x[0,:], x[1,:], '.b', markersize=16)
    plt.legend( h, "grr")
    
    