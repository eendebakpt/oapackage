# -*- coding: utf-8 -*-
"""
Example script to use Nauty from Python

"""

#%% Load packages
import numpy as np
import oapackage; 

# define some graph and colors
G= np.zeros( (5,5), dtype=int); G[0,1]=G[0,2]=G[0,3]=G[1,3]=1
colors = [0,0,0,1,1]

tr = oapackage.reduceGraphNauty(G, colors=colors, verbose=2)
Gx=oapackage.transformGraphMatrix(G, tr)

print('input graph: ')
print(G)

print('reduced graph: ')
print(Gx)


    

    




          