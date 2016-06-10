#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 12:16:43 2013

@author: eendebakpt
"""

#%% Load packages

from __future__ import print_function
import os,sys
import argparse
import numpy as np


try:
    import oapackage
except:
    print('please install the Orthogonal Array package, see https://github.com/eendebakpt/oapackage')
    raise

def filterArray(al, fname, verbose=0):
    if fname=='double' or fname=='dr':
        rs=al.row_symmetry_group()
        if np.any(np.array(rs.gsize)>1):
            return True
        else:
            return False    
    elif fname=='eoc' or fname=='evenoddconf':
        if not oapackage.isConferenceFoldover(al):
            return True
        else:
            return False
    else:
        return True

#%%
#sys.argv += ['test.oa', 'test2.oa']

parser = argparse.ArgumentParser()
parser.add_argument("inputfile")
parser.add_argument("outputfile")
parser.add_argument('-v', '--verbose', default=1, type=int)
parser.add_argument('-f', '--format', help='format for output file' , default='T', type=str)
parser.add_argument('-s', '--filter', help='type of filter to use' , default='double', type=str)
parser.description=('Possible options for the filter: double (select designs with double rows), eoc (even-odd conference matrix)')

args=parser.parse_args()
verbose=args.verbose
fname=args.filter
oaformat=args.format
inputfile=args.inputfile
outputfile=args.outputfile

#%%
if 0:            
    inputfile='/home/eendebakpt/tmp/mee3.oa'
    outputfile='/home/eendebakpt/tmp/test.oa'
    fname='eoc'
    oaformat='T'    
    verbose=2       
#%% Load file

fmode=oapackage.arrayfile_t.parseModeString(oaformat)
infile = oapackage.arrayfile_t(inputfile)
narrays=infile.narrays
#outfile = oapackage.arrayfile_t(outputfile, infile.nrows, infile.ncols, -1, fmode )

if not infile.isopen():
    raise Exception('could not open file %s')

print('oaselect: reading %d arrays from %s' % (infile.narrays, inputfile) )

#%%
outlist=[]

if narrays==-1:
    narrays=oapackage.arrayfile_t.NARRAYS_MAX
for ii in range(narrays):
    oapackage.tprint('oaselect: parse array %d/%d'  % (ii, narrays) )
    al=infile.readnext()
    if al.index==-1:
        break
    if verbose>=3:
        print('  index %d'  % al.index)
    if filterArray(al, fname, verbose=0):
        outlist.append(al)
        #outfile.append_array(al)    
if narrays>0:
    oapackage.tprint('oaselect: parse array %d/%d'  % (ii, narrays), dt=-1 )

oapackage.writearrayfile(outputfile, outlist, fmode, infile.nrows, infile.ncols)

#%%
#outfile.closefile()
#print('oaselect: done (selected %d arrays)' % outfile.narraycounter)
print('oaselect: done (selected %d arrays)' % len(outlist))
