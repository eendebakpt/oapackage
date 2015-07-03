# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 11:47:57 2012

@author: eendebakpt
"""

#%% Load packages
import oapackage
from oapackage.graphtools import designReduceBliss

arrayclass = oapackage.arraydata_t(2, 8, 0, 4)
al=arrayclass.randomarray()

alx, tt = designReduceBliss(al, arrayclass, verbose=1)

alx.showarray()


    

    




          