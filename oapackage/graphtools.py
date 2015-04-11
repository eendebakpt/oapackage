# -*- coding: utf-8 -*-
"""

@author: eendebakpt
"""

#%% Load packages
import sys
import os
import numpy as np
from scipy.special import binom
from math import factorial

#%%
def oa2graph(A, adata, verbose=1):
     """      
     %OA2GRAPH Convert orthogonal array to graph representation
     %
     % The conversion method is as in Ryan and Bulutoglu.
     % The resulting graph is bi-partite.
     % The graph representation can be used for isomorphism testing.
     %    im, colors, r = oa2graph( A, oadata );
    
     """    
     A=np.array(A)
     nrows=adata.N
     ncols=A.shape[1]     
     nRowVertices = A.shape[0]
     nColumnLevelVertices = sum( adata.getS() )
     nVertices = adata.N+ncols+nColumnLevelVertices        
     nColVertices = ncols
     colOffset = adata.N;
    
     s=np.array(adata.getS())
     vertexOffsets = adata.N + ncols + np.hstack( (0, s[0:-1]) ).cumsum() ;
     colors = np.hstack( (np.zeros(adata.N), np.ones(ncols), 2*np.ones(nColumnLevelVertices) ))

     im = np.zeros((nVertices, nVertices)) #  incidence matrix
    
     for row in range(0, nrows):
        idx = A[row, :] + vertexOffsets
        im[row, idx] = 1
        im[idx, row] = 1
    
     if nColVertices>0:
        colidx=2
        for col in range(0, ncols):
            sel = vertexOffsets[col] + range(0, s[col] )
            im[colOffset+col, sel]=colidx
            im[sel, colOffset+col]=colidx

    # The non-row vertices do not have any connections to other non-row vertices.
    
     xy=np.zeros( (2, nrows+s.sum()))

     # calculate positions
     for row in range(0, nrows):
        xy[:, row] = np.array([0,row])

     pos = nrows;
     for col in range(0, ncols):
        for ss in range(0, s[col]):
            xy[:, pos] = np.array([2+ss/s[col], col])
            pos = pos + 1

     return im, colors, dict( {'adata': adata, 'im': im, 'colors': colors, 'nVertices': nVertices})


import oapackage

def graph2arrayTransformation(pp, arrayclass, verbose=0):
    """ From a relabelling of the graph return the corresponding array transformation """
    ppi=np.zeros( len(pp), )
    ppi[pp]=range(len(pp))
    ppi=np.array(ppi).astype(int)
    pp=np.array(pp).astype(int)

    # extract colperms and rowperm and levelperms from this...
    
    rowperm=np.array((pp[0:arrayclass.N]))
    rowperm=rowperm-rowperm.min()
    colperm=np.array( (pp[arrayclass.N:(arrayclass.N+arrayclass.ncols)]) )
    colperm=np.argsort(colperm) # colperm-colperm.min()
    ns=np.sum( arrayclass.getS() )
    lvlperm = np.array( (pp[(arrayclass.N+arrayclass.ncols):(arrayclass.N+arrayclass.ncols+ns)]) )
    lvlperm=lvlperm-lvlperm.min()


    #%%
    ttr=oapackage.array_transformation_t(arrayclass)
    ttr.setrowperm(rowperm)
    ttr=ttr.inverse()
    ttc=oapackage.array_transformation_t(arrayclass)
    ttc.setcolperm(colperm)
    #ttc=ttc.inverse()
    #tt=ttc
    #tt=ttc*ttr
    
    ttl=oapackage.array_transformation_t(arrayclass)
    
    ncols=arrayclass.ncols
    cs=np.hstack( ([0], np.cumsum(arrayclass.getS()) ) )
    lp=[]
    for ii in range(ncols):
        ww=lvlperm[cs[ii]:cs[ii+1]]
        ww=ww-ww.min()
        ww=np.argsort(ww) 
        lp.append( ww)
        ttl.setlevelperm ( ii, ww)
        #ttl.setlevelperm ( colperm[ii], ww)
    
    ttl=ttl.inverse()
    
    tt=ttr*ttc*ttl
    return tt

#

try:
    import PyBliss
except Exception as e:
    print(e)
    print('oapackage.graphtools: could not import PyBliss!')
    pass


def findCanonicalPyBliss(gg, arrayclass, verbose=0):
    """ Find a canonical labelling using PyBliss """
    G = PyBliss.Graph()
    for ii in range(0, arrayclass.N):
        G.add_vertex(ii, color=0)
    
    for ii in range(0, arrayclass.ncols):
        jj=arrayclass.N+ii
        G.add_vertex(jj, color=1)
    
    nn=np.sum(arrayclass.getS())
    for ii in range(0, nn):
        jj=arrayclass.N+arrayclass.ncols+ii
        G.add_vertex(jj, color=2)
    
    if 1:
        ee=zip( list(gg.nonzero()[0]), list(gg.nonzero()[1]) )
        ee=[e for e in ee if e[0]<=e[1]]
        for e in ee:
            G.add_edge(e[0],e[1])
    if verbose:
        print('graph: %d vertices, %d edges' % (G.nof_vertices(), -1 ) )    
    
    canlab = G.canonical_labeling()
    canlab= [canlab[k] for k in sorted(canlab.keys()) ] 
    if verbose:
        print(canlab)
    return canlab
#canlab=findCanonical2(gg, arrayclass, verbose=0)


    
def reduceBliss(al, arrayclass, verbose=1):

    if arrayclass.ismixed():
        print('reduction of mixed classes not implemented?')
        pass
        #raise Exception('reduction of mixed classes not implemented')

    gg, colors, idata2=oa2graph(al, arrayclass )
    
    pp=findCanonicalPyBliss(gg, arrayclass, verbose=0)
    tt=graph2arrayTransformation(pp, arrayclass)

    return pp, tt
    
#%%

def selectIsomorphismClasses(sols, arrayclass, verbose=1):
    """ Select isomorphism classes from a list of designs """
    # perform check on array data type
    mm=[]
    for ii, al in enumerate(sols):
        pp,tt=reduceBliss(al, arrayclass, verbose>=2)
        tt=graph2arrayTransformation(pp, arrayclass)
        alx=tt.apply(al)
        mm.append(np.array(alx))
        pass
        # convert design to graph representation
    
        # convert to canonical form
    
    # perform uniqueness check
    nn=len(mm)
    qq=np.array( [None]*nn, dtype=object)
    for ii in range(nn):
        qq[ii]=mm[ii].flatten()
        
        
    # Trick to make unique work...
    nx=qq[0].size
    dt=qq[0].dtype.descr*nx
    qqq=np.array(qq, dtype=dt)
        
    a,b=np.unique(qqq, return_inverse=True)
    #a,b=np.unique(qqq, return_index=True)
    
    return b, mm

    