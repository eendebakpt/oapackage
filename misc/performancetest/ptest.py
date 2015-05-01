# -*- coding: utf-8 -*-
"""
Python performance test

Log with:
    
    > python ptest.py > test0-ddmmyyy.txt
"""

#%% Load necessary packages 
import sys
import os
#import numpypy
import numpy as np
import time
import glob
import subprocess
import time
import getopt
import platform
from imp import reload

oadir=os.path.join(os.path.expanduser('~'), 'misc/oa/oacode/')
sys.path.append(oadir)
sys.path.append(os.path.join(oadir, 'python'))
   
import oalib; reload(oalib)
import oapackage.oahelper as oahelper; reload(oahelper)
tmp=oalib.log_print(-oalib.SYSTEM, '');

#%% Define test functions

def pytest7(verbose=1):
	""" Reduced array with random transformtion to LMC form """
	al=oalib.exampleArray(3)
	adata=oalib.arraylink2arraydata(al)
	reduction=oalib.LMCreduction_t(adata)

	for ii in range(10):
	    if verbose>=2:
	    	    print('pytest7: randomized reduction of array %d'  % ii)
	    reduction.transformation.randomize()
	    al2=reduction.transformation.apply(al)
	    
	    alr=al2.reduceLMC()
	    c=al==alr
	    if not c:
	       print('pytest7: error: reduction of randomized array failed!' )

#%%

def pytest6(verbose=1):
    """ Check generation of OA(32, t, 2^a) """
    N=32; k=8; t=3; l=[2]*8; rr=[]
    gtrr=[3,5,10,17,33]
    s=oalib.intVector(l)
    adata = oalib.arraydata_t(s, N, t, k)

    if verbose:
        print('pytest6: run different algorithms on the same case' )
    algs=[oalib.MODE_ORIGINAL, oalib.MODE_J4]
    for ii, alg in enumerate(algs):
        algname=oalib.algnames(alg)
        if verbose>=2:
            print('pytest6: running %s, alg %s' % (adata.fullidstr(), algname ))
        rr=[]
        tmp = oahelper.runExtend(N,k,t,l, verbose=verbose, nums=rr, algorithm=alg);
        if not rr==gtrr:
            print('pytest6: case %s' % adata.fullidstr() )
            print('   algorithm %s: error: incorrect number of arrays! %s -> %s' % (algname, rr, gtrr) )

#%%
def pytest3(verbose=2):
    """ Test rank and Defficiency properties """
    k=6
    strength=3
    N=40
    lll=oahelper.runExtend(N,k,strength,2, verbose=0)

    rnk0=np.array([19, 19, 21, 21, 19, 22, 22, 22, 22]) # ground-truth
    D0=np.array([0.0, 0.0, 0.0, 0.1703475035645109, 0.0, 0.8107439985672059, 0.8969911485707478, 0.8953919248315042, 0.9007069559712273]) # gt

    rnk=[ oalib.array2xf(al).rank() for al in lll]
    D=[ al.Defficiency() for al in lll]
    
    if (rnk-rnk0).any():
        if verbose:
            print('  rank test FAILED!')
            return False
    if not (np.abs(D0-D)<1e-6).all():
        if verbose:
            print('  Defficiency test FAILED!')
        return False
    return True
    
def pytest4(verbose=1):
    """ Check generation of OA(16,2, 4^2 2^a) """
    N=16; k=6; t=2; l=[4,4,2]; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    #cases.append('runExtend(N=%d, k=%d, t=%d, l=%s)' %(N,k,t, str(l)) )
    if not rr==[2,7,17,27]:
        print('ERROR: incorrect number of arrays! %s' % rr)

def pytest2level(verbose=1):
    """ Check generation of OA(16,2, 4^2 2^a) """
    N=32; k=10; t=3; l=[2]; rr=[]
    t0=time.time()
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr, algorithm=oalib.MODE_ORIGINAL)
    dt=time.time()-t0
    t0=time.time()
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr, algorithm=oalib.MODE_LMC_2LEVEL)
    dt2level=time.time()-t0

    ll=l*k
    s=oalib.intVector(ll)
 
    adata = oalib.arraydata_t(s, N, t, k)

    if verbose:
        print('case %s: 2-level method %.2f [s] -> %.2f [s]' % (adata.idstr(), dt, dt2level) )
    #cases.append('runExtend(N=%d, k=%d, t=%d, l=%s)' %(N,k,t, str(l)) )


def pytest5(verbose=1):
	""" Performance of analysis algorithms """
	cmd = 'export OMP_NUM_THREADS=1; oaanalyse --gwp -A ../testdata/test64.oa'
	if verbose:
		print('pytest5: running oaanalyse')
	t0=time.time()
	os.system(cmd)
	ta=(time.time()-t0)
	return ta

def pytest2(oadir, verbose=1):
    afile=os.path.join(oadir, '.', 'testdata', 'unstable-rank.oa')
    sols=oalib.readarrayfile(afile)
    afile=os.path.join(oadir, '.', 'testdata', 'unstable-rank-extended.oa')
    sols=oalib.readarrayfile(afile)

    for ii,al in enumerate(sols):
        vv=oalib.doubleVector(3)
        rnk=oalib.array_rank_D_B(al, vv)
    
        al2=oalib.array2xf(al)
        rnk2=al2.rank()
        if (not rnk2==rnk) or (not rnk<29):
            print('warning: problem?!: %d %d %d' %  (ii, rnk, rnk2))
        if verbose>=2:
            print('array %d: rank: %d %d' % (ii, rnk, rnk2))

#    adata=oalib.arraylink2arraydata(al)    
#    oalib.extend_array(al.array, adata, adata.cols)
    
def pytest(verbose=1):
    t0=time.time()
    tt=[]
    cases=[]
    
    N=48; k=6; t=3; l=2; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    tt.append(time.time()-t0)
    cases.append('runExtend(N=%d, k=%d, t=%d, l=%d)' %(N,k,t,l))
    if not rr==[4,10,45]:
        print('ERROR: incorrect number of arrays! %s' % rr)

    N=32; k=6; t=3; l=2; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    tt.append(time.time()-t0)
    cases.append('runExtend(N=%d, k=%d, t=%d, l=%d)' %(N,k,t,l))
    if not rr==[3,5,10]:
        print('ERROR: incorrect number of arrays! %s' % rr)
    
    N=16; k=16; t=2; l=2; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    tt.append(time.time()-t0)
    cases.append('runExtend(N=%d, k=%d, t=%d, l=%d)' %(N,k,t,l))
    if not rr==[3, 5, 11, 27, 55, 80, 87, 78, 58, 36, 18, 10, 5, 0]:
        print('ERROR: incorrect number of arrays! %s' % rr)

    N=96; k=8; t=4; l=2; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    tt.append(time.time()-t0)
    cases.append('runExtend(N=%d, k=%d, t=%d, l=%d)' %(N,k,t,l))
    if not rr==[4, 9, 4, 0]:
        print('ERROR: incorrect number of arrays! %s' % rr)

    dt=time.time()-t0
    if verbose:
        print('total time %.3f [s]' % dt)


    if verbose:
        print('extra testing cases')
    N=24; k=5; t=2; l=[3,2]; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    tt.append(time.time()-t0)
    cases.append('runExtend(N=%d, k=%d, t=%d, l=%s)' %(N,k,t,str(l)))
    if not rr==[4,29,573]:
        print('ERROR: incorrect number of arrays! %s' % rr)

    N=18; k=9; t=2; l=[2,3]; rr=[]
    oahelper.runExtend(N,k,t,l, verbose=1, nums=rr)
    tt.append(time.time()-t0)
    cases.append('runExtend(N=%d, k=%d, t=%d, l=%s)' %(N,k,t, str(l)) )
    if not rr==[3, 15, 48, 19, 12, 3, 0]:
        print('ERROR: incorrect number of arrays! %s' % rr)

    dt2=time.time()-t0
    return (dt, tt, cases, dt2)

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#%%

def testExtendBinary(verbose=1):
    repos=dict()
    repos['oatest1.txt'] = 'result-16.2-2-2-2-2-2-2-2.oa';
    repos['oatest2.txt'] = 'result-18.3-3-3-3-3.oa';
    repos['oatest3.txt'] = 'result-8.4-2-1-4-2-2.oa';
    repos['oatest4.txt'] = 'result-28.7-2-2-2.oa';
    repos['oatest5.txt'] = 'result-16.2-2-2-2-2-2-2-2-2.oa';
    repos['oatest6.txt'] = 'result-25.5-5-5-5.oa';
    repos['oatest7.txt'] = 'result-64.2-2-2-2-2-2.oa';
    repos['oatest8.txt'] = 'result-56.2-2-2-2-2.oa';


#testExtendBinary()

#%%
def main(argv=None):
    """ Main testing function """
    print('OA performance testing')
    ss=oalib.version()
    print('OAlib: version %s' % ss)
    ss=oalib.compile_information()
    print(ss)
    print('System:')
    print('  Python: ' + ' '.join(sys.version.split('\n')))
    print('  Machine: ' + ' '.join(['%s' % x for x in platform.uname()]))

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error as msg:
             raise Usage(msg)
        # more code, unchanged
    except Usage as err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
        
    print('\nRunning tests:\n')
    (dt, tt, cases, dt2)=pytest()    
   
    pytest2(oadir)
    pytest3(verbose=0)
    pytest4(verbose=1)
    ta1 = pytest5(verbose=1)
    pytest6(verbose=0)
    pytest7(verbose=0)

    pytest2level(verbose=1)

    print('-----')
    
    #tt=[];dt=0

    dtx = ta1

    print('\nResults:\n')    
    for ii, t in enumerate(tt):
        print('%s: time %.2f [s]' % (cases[ii], t))
    print('Total time: %.3f [s], %.3f [s], %.3f [s]' % (dt, dt2, dtx))
    print('   should be of order 4.4 [s], 4.6 [s], 5.9 [s] (woelmuis)')
    print('   should be of order 3.2 [s], 3.3 [s], 4.7 [s] (marmot)')
    
if __name__ == "__main__":
    main()    
    #sys.exit(main())
    

    
def timeconfig(configfile = 'oaconfig.txt', timebin='/usr/bin/time', pdir='/home/eendebakpt/misc/oa/oacode/performancetest'):  
    os.chdir(pdir)
    #res,v =os.system('%s --format="%%E %%S %%U" ./oaextendsingle -c %s -l 1' % (timebin, configfile) )
    res = subprocess.check_output([('cd %s;ls;' % pdir) + timebin, '--format="%%E %%S %%U"', './oaextendsingle', '-c %s -l 1' % ( configfile) ])
    return res

    #res = subprocess.check_output([timebin, '--format="usr %%E"', './oaextendsingle', '-l 1'])
    try:
        res =os.system('./oaextendsingle -c %s -l 0 > /dev/null' % (configfile) )
    except:
        print('error!')
    

#res = subprocess.check_output(['time'  , 'ls' ], shell=True)

#import sys; oadir='/home/eendebakpt/misc/oa/oacode/performancetest'; sys.path.append(oadir)
#import copy
#import ctypes
#import ptetools
#from ptetools import *

#from ctypes import *
#oalib = cdll.LoadLibrary("oalib.so")


