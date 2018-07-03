#!/usr/bin/python
"""

Example script for Python interface to Orthogonal Array code.

> pygmentize -f html example_python_testing.py

@author: Pieter Eendebak
"""
#%% Load packages
import sys, os
import numpy as np

print('Python OA Interface Check')
import oalib;
import oapackage.oahelper as oahelper


#%%
print('-- Create design configuration --')
adata=oalib.arraydata_t(2, 32, 2, 7)
oaoptions=oalib.OAextend()
oaoptions.setAlgorithmAuto(adata)
adata.show()

print('-- Create root element --')
al=adata.create_root()
al.show()
al.showarray()

# Extend
print('-- Extend arrays --')

print('Extend to 3 columns')
newsols=oalib.extend_array(al, adata, oaoptions)
print(newsols)

print('Extend to 4 columns')
newsols2=oalib.extend_arraylist(newsols, adata, oaoptions)
print(newsols2)

print('-- Analysing properties of arrays --')
for ii in range(0, 10):
    al=newsols2[ii]
    print('array %d: generalized word length pattern: %s' % (ii, str(al.GWLP()) ) )
    print('          D-efficiency: %.3f' % ((al.Defficiency()) ) )


al=oalib.exampleArray()
al.showarray()
print('D-efficiency %f, rank %d' % (al.Defficiency(), al.rank()) )

#%% Write al array to disk
al=oalib.exampleArray()
r=oalib.writearrayfile('test.oa', al)
oahelper.oainfo('test.oa')

#%% Convert to Numpy array
al=oalib.exampleArray(0)
al.showarray()
al[2,1]
X=al.getarray()
X


