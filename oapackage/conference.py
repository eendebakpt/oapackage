# -*- coding: utf-8 -*-
""" Module with functions to work with conference designs

"""

from __future__ import print_function

import numpy as np
import itertools


import oapackage

# %%


def momentMatrix(k):
    """ Return the moment matrix of a conference design

    Args:
        k (int): number of columns of the conference design
    Returns:
        array: moment matrix
    """
    pk = int(1 + 0.5 * k * (k + 1) + k)
    M = np.zeros((pk, pk))
    M[0, 0] = 1
    M[0, int(0.5 * k * (k + 1) + 1):] = 1. / 3
    M[int(0.5 * k * (k + 1) + 1):, 0] = 1. / 3
    M[1:(k + 1), 1:(k + 1)] = np.eye(k) * 1. / 3
    M[(k + 1):int(0.5 * k * (k + 1) + 1), (k + 1):int(0.5 * k * (k + 1) + 1)] = np.eye(int(0.5 * k * (k - 1))) / 9.
    M[int(0.5 * k * (k + 1) + 1):, int(0.5 * k * (k + 1) + 1):] = np.eye(k) / 5 + (np.ones((k, k)) - np.eye(k)) / 9.
    return M


def _leftDivide(A, B):
    r""" Perform left division of a matrix

    Args:
        A (aray)
        B (array)
    Returns:
        array: the result of A\B
    """
    solution = np.linalg.solve(A.T.dot(A), A.T.dot(B))
    return solution


def modelStatistics(dsd, verbose=0, moment_matrix=None, use_condition_number=True):
    """ Calculate statistics of a definitive screening design from the model matrix

    Args:
        dsd (array): definitive screening design
    Returns:
      list: calculated statistics. Calculated are whether the model is estible, the Defficiency and the inverse APV.
    """
    ncolumns = dsd.shape[1]
    modelmatrix = np.array(oapackage.conference_design2modelmatrix(dsd, 'q', verbose=0))
    M = (modelmatrix.T).dot(modelmatrix)

    if use_condition_number:
        fullrank = np.linalg.cond(M) < 1000
    else:
        mr = np.linalg.matrix_rank(M)
        fullrank = (mr == modelmatrix.shape[1])
        if verbose >= 2:
            print('modelStatistics: fullrank %d, condition number %.4f' % (fullrank, np.linalg.cond(M), ))

    if verbose >= 2:
        print('modelStatistics: condition number: %s' % (np.linalg.cond(M)))
    if fullrank:
        if moment_matrix is None:
            moment_matrix = momentMatrix(ncolumns)
        Eest = 1
        pk = int(1 + ncolumns + ncolumns * (ncolumns + 1) / 2)
        kappa = np.linalg.det(M)
        lnkappa = np.log(kappa) / (pk)
        Defficiency = np.exp(lnkappa)

        apv = np.trace(_leftDivide(M, moment_matrix))
        invAPV = 1 / apv
    else:
        Eest = 0
        Defficiency = 0
        invAPV = 0

    return Eest, Defficiency, invAPV


def conferenceProjectionStatistics(array, ncolumns=4, verbose=0):
    """ Calculate the projection statistics of a conference design

    The PECk, PICk and PPCk are calculated with k the number of columns specified.
    The projection statistics are calculated by determined the :meth:`modelStatistics` of all k-column subdesigns
    and then taking the mean of the statistics for the subdesigns.
    For more details of the calculation, see https://oapackage.readthedocs.io/.

    Args:
        array (array): conference design
        ncolumns (int): number of columns on which to project

    Returns:
        pec, pic, ppc (float): calculated statistics
    """
    nc = array.shape[1]

    number_combinations = oapackage.choose(nc, ncolumns)
    Eestx = np.zeros(number_combinations)
    Deff = np.zeros(number_combinations)
    invAPV_values = np.zeros(number_combinations)
    dsd = oapackage.conference2DSD(oapackage.array_link(array))
    moment_matrix = momentMatrix(ncolumns)
    for idx, c in enumerate(list(itertools.combinations(range(nc), ncolumns))):
        proj_dsd = dsd.selectColumns(c)
        Eest, D, invAPV = modelStatistics(proj_dsd, verbose=0, moment_matrix=moment_matrix)

        Deff[idx] = D
        Eestx[idx] = Eest
        invAPV_values[idx] = invAPV
    pec, pic, ppc = np.mean(Eestx), np.mean(Deff), np.mean(invAPV_values)
    if verbose:
        print('conferenceProjectionStatistics: projection to %d columns: PEC %.3f PIC %.3f PPC %.3f  ' %
              (ncolumns, pec, pic, ppc))
    return pec, pic, ppc
