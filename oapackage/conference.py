# -*- coding: utf-8 -*-
""" Module to generate D-optimal designs

For more information see: https://doi.org/10.1080/00401706.2016.1142903

Pieter Eendebak <pieter.eendebak@gmail.com>

"""

from __future__ import print_function

import os
import numpy as np
import time
import itertools


import oapackage
import oapackage.markup as markup
import oapackage.oahelper as oahelper


#%%

def momentMatrix(k):
    """ Return the moment matrix of a conference design 

    Args
        k (int): number of columns of the conference design 
    Returns
        array: moment matrix
    """
    pk = int(1 + 0.5 * k * (k + 1) + k)
    M = np.zeros((pk, pk))
    M[0, 0] = 1
    M[0, int(0.5 * k * (k + 1) + 1):] = 1 / 3
    M[int(0.5 * k * (k + 1) + 1):, 0] = 1 / 3
    M[1:(k + 1), 1:(k + 1)] = np.eye(k) * 1. / 3
    M[(k + 1):int(0.5 * k * (k + 1) + 1), (k + 1):int(0.5 * k * (k + 1) + 1)] = np.eye(int(0.5 * k * (k - 1))) / 9
    M[int(0.5 * k * (k + 1) + 1):, int(0.5 * k * (k + 1) + 1):] = np.eye(k) / 5 + (np.ones((k, k)) - np.eye(k)) / 9
    return M


def _leftDivide(A, B):
    """ Perform left division of a matrix

    Args:
        A (aray)
        B (array)
    Returns:
        array: the result of A\B
    """
    x, resid, rank, s = np.linalg.lstsq(A, B, rcond=None)
    #x = lin.solve(A.T.dot(A), A.T.dot(B))
    return x


def modelStatistics(dsd, verbose=0):
    """ Calculate statistics of a definitive screening design from the model matrix

    Args:
        dsd (array): definitive screening design
    """
    nrows = dsd.shape[0]
    ncolumns = dsd.shape[1]
    modelmatrix = oapackage.array2modelmatrix(dsd, 'q')
    modelmatrix = oapackage.array_link(modelmatrix.astype(int))

    Eest = 0
    M = (np.array(modelmatrix).T).dot(np.array(modelmatrix))
    mr = np.linalg.matrix_rank(M)

    if verbose >= 2:
        print('conferenceProjectionStatistics: condition number: %s' % (np.linalg.cond(M)))
    if verbose:
        print('%d, cond %.5f' % (mr == modelmatrix.shape[1], np.linalg.cond(M), ))
    if mr == modelmatrix.shape[1]:  # np.linalg.cond(np.array(A).dot(np.array(A).T))<1000:
        Eest = 1
        pk = int(1 + ncolumns + ncolumns * (ncolumns + 1) / 2)
        kappa = np.linalg.det(M)
        lnkappa = np.log(kappa) / (pk)
        Defficiency = np.exp(lnkappa)

        apv = np.trace(_leftDivide(M, momentMatrix(ncolumns)))
        invAPV = 1 / apv
    else:
        Eest = 0
        Defficiency = 0
        invAPV = 0

    return Eest, Defficiency, invAPV


def conferenceProjectionStatistics(al, ncolumns=4, verbose=0):
    """ Calculate the projection statistics of a conference design

    The PECk, PICk and PPCk are calculated with k the number of columns specified.

    Args:
        al (array): conference design
        ncolumns (int): number of column on which to project

    Returns:
        pec, pic, ppc (float)
    """
    nc = al.shape[1]
    AA = np.array(al)
    Eestx = []
    Deff = []
    invAPV_values = []
    for c in list(itertools.combinations(range(nc), ncolumns)):
        X = AA[:, c]
        dsd = oapackage.conference2DSD(oapackage.array_link(X))
        k = X.shape[1]

        Eest, D, invAPV = modelStatistics(dsd, verbose=0)

        Deff += [D]
        Eestx += [Eest]
        invAPV_values += [invAPV]
    pec, pic, ppc = np.mean(Eestx), np.mean(Deff), np.mean(invAPV_values)
    if verbose:
        print('conferenceProjectionStatistics: projection to %d columns: PEC %.3f PIC %.3f PPC %.3f  ' % (ncolumns, pec, pic, ppc))
    return pec, pic, ppc

