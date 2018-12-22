# -*- coding: utf-8 -*-
""" Module to work with designs and graphs

@author: Pieter Eendebak <pieter.eendebak@gmail.com>
"""

# %% Load packages
import numpy as np
import oapackage

# %%


def oa2graph(al, adata, verbose=1):
    """
    Convert orthogonal array to graph representation

    The conversion method is as in Ryan and Bulutoglu. The resulting graph is bi-partite.
    The graph representation can be used for isomorphism testing.

    Example
    -------
    >     im, colors, r = oa2graph( A, oadata );

    """
    A = al.getarray(verbose=0)
    if verbose:
        print('oa2graph: array of shape %s' % (A.shape, ))
    nrows = adata.N
    ncols = A.shape[1]
    nColumnLevelVertices = sum(adata.factor_levels())
    nVertices = adata.N + ncols + nColumnLevelVertices
    nColVertices = ncols
    colOffset = adata.N

    s = np.array(adata.factor_levels())
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
    return im, colors, dict({'adata': adata, 'im': im, 'colors': colors, 'nVertices': nVertices})


def selectIsomorphismClasses(sols, verbose=1):
    """ Select isomorphism classes from a list of designs

    Args:
        sols (list of arrays): list of arrays from which to determine the unique ones
        verbose (int): verbosity level
    Return:
        indices (list of integers): indices of the isomorphism classes
        mm (list of arrays): the arrays in normal form


    Example:
        >>> import oapackage.graphtools; import numpy as np
        >>> sols=[oapackage.exampleArray(idx) for idx in [26,26,27,26]]
        >>> idx, mm = oapackage.graphtools.selectIsomorphismClasses(sols, verbose=0)
        >>> _, unique_indices = np.unique(idx, return_index=True)
        >>> print('found %d isomorphism classes from %d arrays' % (len(unique_indices), len(sols)) )
        found 2 isomorphism classes from 4 arrays
        >>> representatives = [sols[idx] for idx in unique_indices]

    """

    mm = []
    for ii, al in enumerate(sols):
        if verbose:
            oapackage.tprint(
                'selectIsomorphismClasses: process aray %d/%d' % (ii, len(sols)))
        al = oapackage.makearraylink(al)

        tt = oapackage.reduceOAnauty(al, verbose >= 2)

        alx = tt.apply(al)
        mm.append(np.array(alx))

    # perform uniqueness check
    nn = len(mm)
    qq = np.array([None] * nn, dtype=object)
    for ii in range(nn):
        qq[ii] = mm[ii].flatten()

    # Trick to make unique work...
    _, indices = np.unique(np.vstack(qq), axis=0, return_inverse=True)

    if verbose >= 1:
        print('selectIsomorphismClasses: reduce %d to %d' %
              (len(sols), np.unique(indices).size))

    return indices, mm


def test_select_isomorphism():
    ll = [oapackage.exampleArray(0), oapackage.exampleArray(0)]
    indices, mm = selectIsomorphismClasses(ll, verbose=3)
    assert(indices[0] == indices[1])
    assert(len(mm) == 2)
    assert(np.all(mm[0] == mm[1]))
