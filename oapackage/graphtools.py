# -*- coding: utf-8 -*-
"""

@author: eendebakpt
"""

#%% Load packages
import numpy as np
import oalib
import oapackage
import sys


#%%


def oa2graph(al, adata, verbose=1):
    """
    Convert orthogonal array to graph representation

    The conversion method is as in Ryan and Bulutoglu.
    The resulting graph is bi-partite.
    The graph representation can be used for isomorphism testing.

    Example
    -------
    >     im, colors, r = oa2graph( A, oadata );

    """
    #A = np.array(al)
    A = al.getarray(verbose=0)
    #print('type A: %s' % type(A) )
    nrows = adata.N
    ncols = A.shape[1]
    nRowVertices = A.shape[0]
    nColumnLevelVertices = sum(adata.getS())
    nVertices = adata.N + ncols + nColumnLevelVertices
    nColVertices = ncols
    colOffset = adata.N

    s = np.array(adata.getS())
    sc = np.cumsum(s)
    sc0 = np.hstack(([0], sc))
    qq = np.ones(nColumnLevelVertices)  # colors for column level vertices
    for ii in range(sc.size):
        qq[sc0[ii]:sc0[ii + 1]] = 2 + ii
    qq = 2 * np.ones(nColumnLevelVertices)  # colors for column level vertices

    vertexOffsets = adata.N + ncols + np.hstack((0, s[0:-1])).cumsum()
    colors = np.hstack(
        (np.zeros(adata.N), np.ones(ncols), qq))

    im = np.zeros((nVertices, nVertices))  # incidence matrix

    for row in range(0, nrows):
        idx = A[row, :] + vertexOffsets
        im[row, idx] = 1
        im[idx, row] = 1

    if nColVertices > 0:
        colidx = 2
        for col in range(0, ncols):
            sel = vertexOffsets[col] + range(0, s[col])
            im[colOffset + col, sel] = colidx
            im[sel, colOffset + col] = colidx

    # The non-row vertices do not have any connections to other non-row
    # vertices.

    if 0:
        # xy is an old data structure which is not used any moref
        xy = np.zeros((2, nrows + s.sum()))

        # calculate positions
        for row in range(0, nrows):
            xy[:, row] = np.array([0, row])

        pos = nrows
        for col in range(0, ncols):
            for ss in range(0, s[col]):
                xy[:, pos] = np.array([2 + ss / s[col], col])
                pos = pos + 1

    return im, colors, dict({'adata': adata, 'im': im, 'colors': colors, 'nVertices': nVertices})


import oapackage


def graph2arrayTransformation(pp, arrayclass, verbose=0):
    """ From a relabelling of the graph return the corresponding array transformation """
    ppi = np.zeros(len(pp), )
    ppi[pp] = range(len(pp))
    ppi = np.array(ppi).astype(int)
    pp = np.array(pp).astype(int)

    # extract colperms and rowperm and levelperms from this...

    rowperm = np.array((pp[0:arrayclass.N]))
    rowperm = rowperm - rowperm.min()
    colperm = np.array((pp[arrayclass.N:(arrayclass.N + arrayclass.ncols)]))
    colperm = np.argsort(colperm)  # colperm-colperm.min()
    ns = np.sum(arrayclass.getS())
    lvlperm = np.array(
        (pp[(arrayclass.N + arrayclass.ncols):(arrayclass.N + arrayclass.ncols + ns)]))
    lvlperm = lvlperm - lvlperm.min()

    #%%
    ttr = oapackage.array_transformation_t(arrayclass)
    ttr.setrowperm(rowperm)
    ttr = ttr.inverse()
    ttc = oapackage.array_transformation_t(arrayclass)
    ttc.setcolperm(colperm)

    ttl = oapackage.array_transformation_t(arrayclass)

    ncols = arrayclass.ncols
    cs = np.hstack(([0], np.cumsum(arrayclass.getS())))
    lp = []
    for ii in range(ncols):
        ww = lvlperm[cs[ii]:cs[ii + 1]]
        ww = ww - ww.min()
        ww = np.argsort(ww)
        lp.append(ww)
        ttl.setlevelperm(ii, ww)

    ttl = ttl.inverse()

    tt = ttr * ttc * ttl
    return tt

#


#%%

from oapackage import makearraylink


def selectIsomorphismClasses(sols, verbose=1):
    """ Select isomorphism classes from a list of designs 
    
    Args:
        sols (list of arrays)
        verbose (int)
    Return:
        indices (list)
        
    """

    # perform check on array data type
    mm = []
    for ii, al in enumerate(sols):
        if verbose:
            oapackage.tprint('selectIsomorphismClasses: process aray %d/%d'  % ( ii, len(sols)) )
        al = oapackage.makearraylink(al)

        tt = oapackage.reduceOAnauty(al, verbose >= 2)

        #pp, tt = reduceBliss(al, arrayclass, verbose >= 2)
        #tt = graph2arrayTransformation(pp, arrayclass)
        alx = tt.apply(al)
        mm.append(np.array(alx))
        pass

    # perform uniqueness check
    nn = len(mm)
    qq = np.array([None] * nn, dtype=object)
    for ii in range(nn):
        qq[ii] = mm[ii].flatten()

    # Trick to make unique work...
    nx = qq[0].size
    dt = qq[0].dtype.descr * nx
    qqq = np.array(qq, dtype=dt)

    a, indices = np.unique(qqq, return_inverse=True)

    if verbose >= 1:
        print('selectIsomorphismClasses: reduce %d to %d' %
              (len(sols), np.unique(indices).size))

    return indices, mm


def test_select_isomorphism():
    ll = [oapackage.exampleArray(0), oapackage.exampleArray(0)]
    selectIsomorphismClasses(ll)
