# -*- coding: utf-8 -*-
""" Module to generate D-optimal designs

For more information see: https://doi.org/10.1080/00401706.2016.1142903

Pieter Eendebak <pieter.eendebak@gmail.com>

"""

from __future__ import print_function

import os
import numpy as np
import time
import logging
import warnings

import oalib

import oapackage.markup as markup
import oapackage.oahelper as oahelper
from oapackage.markup import oneliner as e


try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm
except BaseException:
    matplotlib = None


class MissingMatplotLibException(Exception):
    pass


# %%


def array2Dtable(array_list, verbose=1, titlestr=None):
    """ Generate HTML table with information about for a list of designs

    Args:
        array_list (list): list of arrays
        verbose (int): verbosity level
        titlestr (str): Not used
    """
    page = markup.page()
    page.table(style=' border-collapse: collapse;')
    page.tr(style='font-weight: bold; border-bottom: solid 1px black;')
    page.th('Array', style='padding-right:30px; ')
    page.th(('D-efficiency', 'Ds-efficiency', 'D1-efficiency'),
            style='padding-right:14px;')
    page.th(('GWLP'), style='padding-right:14px;')
    page.tr.close()
    for ii, al in enumerate(array_list):
        aidx = ii
        (D, Ds, D1) = al.Defficiencies()
        gwlp = al.GWLP()
        page.tr(style='font-weight: normal;')
        page.td('%d' % aidx, style='padding-right:10px;')
        for statistic in [D, Ds, D1]:
            page.td('%.4f' % statistic, style='padding-right:1px;')
        gstr = oahelper.gwlp2str(gwlp)
        page.td(e.small(gstr), style='padding-right:1px;')
        page.tr.close()
    page.table.close()
    return page

# %%


def generateDscatter(dds, second_index=0, first_index=1, lbls=None, ndata=3, nofig=False, fig=20,
                     scatterarea=80, verbose=0, setWindowRectangle=False):
    """ Generate scatter plot for D and Ds efficiencies

    Args:
        dds (array): array with D-efficiencies
    Returns:
        dict: contains handles to plotting elements
    """
    if matplotlib is None:
        raise MissingMatplotLibException
    data = dds.T
    pp = oahelper.createPareto(dds)
    paretoidx = np.array(pp.allindices())

    nn = dds.shape[0]
    area = scatterarea * np.ones(nn,) / 2
    area[np.array(pp.allindices())] = scatterarea
    alpha = 1.0

    if dds.shape[1] > ndata:
        colors = dds[:, ndata]
    else:
        colors = np.zeros((nn, 1))

    idx = np.unique(colors).astype(int)

    if verbose:
        print('generateDscatter: unique colors: %s' % (idx, ))
    ncolors = idx.size
    try:
        import brewer2mpl
        ncolors = max(ncolors, 4)
        mycmap = brewer2mpl.get_map('Set1', 'qualitative', ncolors).mpl_colors
    except BaseException:
        mycmap = [matplotlib.cm.jet(ii) for ii in np.linspace(0, 1, ncolors)]

    nonparetoidx = np.setdiff1d(range(nn), paretoidx)
    if lbls is None:
        lbls = ['%d' % i for i in range(len(idx))]

    if fig is not None:
        figh = plt.figure(fig)
        plt.clf()
        figh.set_facecolor('w')
        ax = plt.subplot(111)

        ax.scatter(data[first_index, nonparetoidx], data[second_index, nonparetoidx], s=.33 * scatterarea,
                   c=(.5, .5, .5), linewidths=0, alpha=alpha, label='Non-pareto design')

        for jj, ii in enumerate(idx):
            gidx = (colors == ii).nonzero()[0]
            gp = np.intersect1d(paretoidx, gidx)

            color = mycmap[jj]
            cc = [color] * len(gp)
            if verbose:
                print('index %d: %d points' % (ii, gidx.size))
            ax.scatter(data[first_index, gp], data[second_index, gp], s=scatterarea, c=cc,
                       linewidths=0, alpha=alpha, label=lbls[jj])
            plt.draw()

        if data[second_index, :].std() < 1e-3:
            y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            ax.yaxis.set_major_formatter(y_formatter)

        xlabelhandle = plt.xlabel('$D_s$-efficiency', fontsize=16)
        plt.ylabel('D-efficiency', fontsize=16)

        if setWindowRectangle:
            try:
                oahelper.setWindowRectangle(10, 10, 860, 600)
            except Exception as ex:
                print('generateDscatter: setWindowRectangle failed')
                logging.exception(ex)

        plt.axis('image')
        pltlegend = ax.legend(loc=3, scatterpoints=1)  # , fontcolor=almost_black)
        if not nofig:
            plt.show()
        ax.grid(b=True, which='both', color='0.85', linestyle='-')
        ax.set_axisbelow(True)

        if nofig:
            plt.close(figh.number)
        else:
            plt.draw()
            plt.pause(1e-3)
    else:
        ax = None
        xlabelhandle = None
        pltlegend = None
    hh = dict({'ax': ax, 'xlabelhandle': xlabelhandle, 'pltlegend': pltlegend})
    return hh


def generateDpage(outputdir, arrayclass, dds, allarrays, fig=20, optimfunc=(1, 0, 0),
                  nofig=False, urlprefix='', makeheader=True, verbose=1, lbls=None):
    """ Helper function to generate web page with D-optimal design results """
    if verbose:
        print('generateDpage: dds %s' % str(dds.shape))

    pp = oahelper.createPareto(dds)

    narrays = dds.shape[0]
    npareto = pp.number()

    if verbose:
        print('generateDpage: narrays %d' % narrays)

    xstr = oahelper.series2htmlstr(arrayclass, case=1)
    xstrplain = oahelper.series2htmlstr(arrayclass, html=0, case=1)

    if verbose:
        print('generateDpage: selectParetoArrays ')
    paretoarrays = oahelper.selectParetoArrays(allarrays, pp)
    at = array2Dtable(paretoarrays, verbose=1)

    if verbose:
        print('generateDpage: write file with Pareto arrays')

    pfile0 = 'paretoarrays.oa'
    pfile = os.path.join(outputdir, pfile0)
    oalib.writearrayfile(pfile, paretoarrays)

    istrlnk = markup.oneliner.a('paretoarrays.oa', href=urlprefix + pfile0)

    if lbls is None:
        lbls = ['Optimization of $D$']
    if fig is not None:
        hh = generateDscatter(dds, lbls=lbls, fig=fig, nofig=nofig)
        oahelper.niceplot(hh.get('ax', None), despine=True, legend=hh['pltlegend'])

        scatterfile = os.path.join(outputdir, 'scatterplot.png')
        if verbose:
            print('generateDpage: writen scatterplot to %s' % scatterfile)
        plt.savefig(scatterfile, bbox_inches='tight', pad_inches=0.25, dpi=160)

    page = markup.page()

    if makeheader:
        page.init(title="Class %s" % xstrplain,
                  css=('../oastyle.css'),
                  lang='en', htmlattrs=dict({'xmlns': 'http://www.w3.org/1999/xhtml', 'xml:lang': 'en'}),
                  header="<!-- Start of page -->",
                  bodyattrs=dict({'style': 'padding-left: 3px;'}),
                       doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
                       metainfo=({'text/html': 'charset=utf-8', 'keywords': 'orthogonal arrays designs',
                                  'robots': 'index, follow', 'description': 'Even-Odd arrays'}),
                       footer="<!-- End of page -->")

    page.h1('Results for array class %s ' % xstr)

    ss = r'The Pareto optimaly was calculated according to the statistics D, D<sub>1</sub> and D<sub>s</sub>.'
    if npareto == 1:
        page.p('Generated %d arrays, %d is Pareto optimal. %s' %
               (narrays, npareto, ss))
    else:
        page.p('Generated %d arrays, %d are Pareto optimal. %s' %
               (narrays, npareto, ss))

    if narrays > 0:

        scores = calcScore(dds, optimfunc)
        _, dd, sols = selectDn(scores, dds, allarrays, nout=1)
        A = sols[0]
        bestdesignfile = os.path.join(outputdir, 'best-design.oa')
        oalib.writearrayfile(bestdesignfile, A)

        page.h2('Best design')
        page.p('The best design: %s.' %
               e.a('best-design.oa', href=os.path.join(urlprefix, 'best-design.oa')))

        page.p()

        dd = dd[0]
        page.span('D-efficiency: %.4f, ' % A.Defficiency())  # page.br()
        page.span('D<sub>s</sub>-efficiency: %.4f, ' % dd[1])  # page.br()
        page.span('D<sub>1</sub>-efficiency: %.4f' % dd[2])
        page.br()
        page.span('A-efficiency: %.3f' % A.Aefficiency())
        page.br()
        gwlp = A.GWLP()
        # gwlp=','.join(['%.3f' % xq for xq in gwlp])
        gwlp = oahelper.gwlp2str(gwlp, jstr=', ')
        page.span('Generalized wordlength pattern: %s' % gwlp)
        page.br()
        # page.p('D-efficiency: %.3f' % A.Defficiency() )
        pec = oalib.PECsequence(A)
        pec = ','.join(['%.3f' % xq for xq in pec])
        page.span('PEC-sequence: %s' % pec)
        page.br()

    page.h2('Table of Pareto optimal arrays ')

    page.span(str(at))
    page.p('All Pareto optimal arrays: %s' % istrlnk)

    page.img(src=urlprefix + 'scatterplot.png',
             style="margin: 10px; width:95%; min-width: 300px;  max-width:1100px; height: auto; ")

    citationstr = markup.oneliner.a(
        'Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays', href='https://doi.org/10.1002/jcd.20236')

    page.br(clear='both')
    page.p(
        'Citation notice: if you make use of the results on this page, please cite the following paper:')
    page.p('%s, Journal of Combinatorial Designs, Volume 18, Issue 2, pages 123-140, 2010.' %
           citationstr)

    page.p('Generated with oapackage %s, date %s.' %
           (oalib.version(), oahelper.timeString()))

    outfile = os.path.join(outputdir, 'Dresults.html')
    fid = open(outfile, 'wt')
    fid.write(str(page))
    fid.close()
    print('written to file %s' % outfile)

    return outfile


def scoreDn(dds, optimfunc):
    """ Calculate scores from various efficiencies

    Args:
        dds (array): calculated D-efficiencies
        optimfunc (list): parameters for optimization function
    Returns:
        scores (array)
    """
    scores = np.array([oalib.scoreD(dd, optimfunc) for dd in dds])
    return scores


def calcScore(dds, optimfunc):
    """ Calculate D-efficiency score using multiple efficiencies and a weight factor

    Args:
        dds (array): the rows contains the effieciencies for various designs
        optimfunc (array): aray with the weight factors for the efficiencies
    """
    scores = np.array(dds).dot(np.array(optimfunc))

    return scores


def optimDeffPython(A0, arrayclass=None, niter=10000, nabort=2500, verbose=1, alpha=[1, 0, 0], method=0):
    """ Optimize array using specified optimization method

    Args:
        A0 (array_link): design to optimize
        arrayclass (object): contains class of designs to optimize
        alpha (list): specifies the optimization function

    Returns:
        d (array): efficiencies
        A (array): optimized design
    """
    # get factor levels
    if arrayclass is None:
        s = A0.getarray().max(axis=0) + 1
    else:
        s = arrayclass.factor_levels()
    sx = tuple(s.astype(np.int64))
    sx = tuple([int(v) for v in sx])

    sg = oalib.symmetry_group(sx)
    gidx = tuple(sg.gidx)
    gsize = tuple(sg.gsize)
    gstart = tuple(sg.gstart)

    nx = 0

    Dinitial = A0.Defficiency()
    if verbose:
        print('optimDeff: initial Deff %.4f' % Dinitial)
    N = A0.n_rows
    k = A0.n_columns

    alpha_is_function = str(type(alpha)) == "<type 'function'>" or callable(alpha)
    # initialize score
    if alpha_is_function:
        efficiencies = alpha(A0)
    else:
        D, Ds, D1 = A0.Defficiencies()
        efficiencies = alpha[0] * D + alpha[1] * Ds + alpha[2] * D1
    A = A0.clone()
    lc = 0
    for ii in range(0, niter):
        r = np.random.randint(N)
        c = np.random.randint(k)
        r2 = np.random.randint(N)
        # make sure c2 is in proper column
        c2 = gstart[gidx[c]] + oalib.fastrand() % gsize[gidx[c]]

        o = A._at(r, c)
        o2 = A._at(r2, c2)  # no extra error checking
        # swap
        if method == oalib.DOPTIM_SWAP:
            # swap
            if o == o2:
                continue
            A._setvalue(r, c, o2)
            A._setvalue(r2, c2, o)
        elif method == oalib.DOPTIM_NONE:
            pass
        else:
            # flip
            A._setvalue(r, c, 1 - o)

        if alpha_is_function:
            dn = alpha(A)
        else:
            D, Ds, D1 = A.Defficiencies()
            nx = nx + 1
            dn = alpha[0] * D + alpha[1] * Ds + alpha[2] * D1
        if dn >= efficiencies:
            if dn > efficiencies:
                lc = ii
                efficiencies = dn
            if verbose >= 2:
                print('ii %d: %.6f -> %.6f' % (ii, efficiencies, dn))

        else:
            # restore to original
            if method == oalib.DOPTIM_SWAP:
                A._setvalue(r, c, o)
                A._setvalue(r2, c2, o2)
            elif method == oalib.DOPTIM_NONE:
                pass
            else:
                A._setvalue(r, c, o)
        if (ii - lc) > nabort:
            if verbose:
                print('optimDeff: early abort ii %d, lc %d' % (ii, lc))
            break

    if verbose:
        Dfinal = A.Defficiency()
        if Dfinal > Dinitial:
            print('optimDeff: final Deff improved: %.4f -> %.4f' %
                  (Dinitial, A.Defficiency()))
        else:
            print('optimDeff: final Deff %.4f' % A.Defficiency())

    return efficiencies, A


# %%
def filterPareto(scores, dds, arrays, verbose=0):
    """ From a list of designs select only the pareto optimal designs

    Args:
        scores (array): array of scores
        dds (array): array with D-efficiency values
        arrays (list): list of designs
    Returns:
        pareto_scores (list) : list of scores of pareto optimal designs
        pareto_efficiencies (list): list of D-efficiency values of pareto optimal designs
        pareto_designs (list) : list of selected designs
    """
    pp = oahelper.createPareto(dds, verbose=0)
    paretoidx = np.array(pp.allindices())

    pareto_scores = scores[paretoidx]
    pareto_efficiencies = dds[paretoidx]
    pareto_designs = [arrays[i] for i in paretoidx]

    return pareto_scores, pareto_efficiencies, pareto_designs


def selectDn(scores, dds, sols, nout=1, sortfull=True):
    """ Select best arrays according to given scores

    The resulting data is sorted

    Parameters
    ----------
    scores : array
        scores for the designs
    dds : array
        array with efficiencies of the designs
    sols : list
        list of designs
    nout : integer or None
        Number of results to return. None means return all results
    sortfull : boolean
        If True, then sort on both the scores and the dds values

    Returns
    -------
    scores, dds, sols : sorted arrays
    """

    if len(sols) == 0:
        return scores, dds, sols

    if sortfull:
        full_sorting_data = [dds[:, ii] for ii in list(range(dds.shape[1]))[::-1]] + [np.array(scores).flatten()]
        idx = np.lexsort(full_sorting_data)[::-1]
        idx = np.lexsort(np.array(full_sorting_data))[::-1]
    else:
        idx = np.argsort(-scores.flatten())
    scores = scores[idx]
    dds = dds[idx, :]
    sols = [sols[ii] for ii in idx]
    if nout is not None:
        # sort the arrays
        nout = np.minimum(nout, scores.size)
        scores = scores[0:nout]
        dds = dds[range(nout), :]
        sols = sols[0:nout]
    return scores, dds, sols


def Doptimize(arrayclass, nrestarts=10, optimfunc=[
              1, 0, 0], verbose=1, maxtime=180, selectpareto=True, nout=None, method=oalib.DOPTIM_UPDATE, niter=100000, nabort=0, dverbose=1):
    """ Calculate D-efficient designs

    The method uses a coordinate-exchange algorithm find a D-efficient (sometimes called D-optimal) design in the class specified by the
    arrayclass. The optimality is defined in terms of the optimization parameters. The optimization is performed
    multiple times (specified by the nrestarts parameter) to prevent finding a design in a local minmum of the
    target function.


    Args:
      arrayclass (object): Specifies the type of design to optimize
      nrestarts (int): Number of restarts of the algorithm
      optimfunc (list with 3 floats): Gives the optimization weights :math:`\\alpha` of the target function :math:`\\alpha[0] D+\\alpha[1] D_s+\\alpha[2] D_1`
      verbose (int): Verbosity level. A higher numer gives more output
      maxtime (float): Maximum running time of the algorithm. If this time is exceeded the algorithm is aborted.
      selectpareto (bool): default is True. If True then only the Pareto optimal designs are returned
      nout (int or None): Maximum number of designs to return. If None,  return all designs
      method (coordinate_exchange_method_t): Specifies the method use for updating elements in the coordinate-exchange algorithm.
      niter (int): Maximum number of iterations of the coordinate-exchange algorithm
      nabort (int): If no improvements have been found after this number of updates, then abort this run
      dverbose (int): Verbosity level passed to the C++ Doptimize function

    Returns
    -------
            scores : list
                list of scores
            dds: array
                array with calculated efficiencies
            designs: list
                list of generated designs
            nrestarts: int
                number of restarts used


    The optimization target and the Pareto optimality are defined in terms of the D-efficiency, main effect robustness
    (or Ds-optimality) and the D1-efficiency of the design. A full definition of these efficiencies is available
    in the documentation at https://oapackage.readthedocs.io/en/latest/properties.html#optimality-criteria-for-d-efficient-designs. For more details and motivation of these efficiencies,
    see the paper "Two-Level Designs to Estimate All Main Effects and Two-Factor Interactions", https://doi.org/10.1080/00401706.2016.1142903.


    """
    if arrayclass.strength != 0:
        warnings.warn('Doptimize can only handle designs with strength 0', UserWarning)

    if verbose:
        print('Doptim: optimization class %s' % arrayclass.idstr())
    t0 = time.time()

    if optimfunc is None:
        optimfunc = [1, 2, 0]

    if isinstance(optimfunc, list):
        rr = oalib.Doptimize(arrayclass, nrestarts, alpha=optimfunc, verbose=dverbose,
                             method=method, niter=niter, maxtime=maxtime, nabort=nabort)
        dds, sols = rr.dds, rr.designs
        dds = np.array([x for x in dds])
        # needed because of SWIG wrapping of struct type
        sols = [design.clone() for design in sols]
        nrestarts = rr.nrestarts
        scores = np.array(
            [oalib.scoreD(A.Defficiencies(), optimfunc) for A in sols])

        if verbose >= 3:
            print('Doptimize: max score %.3f, max D: %.6f' %
                  (np.max(scores), np.max([A.Defficiency() for A in sols])))
    else:
        # assume optimfunc is a function
        scores = np.zeros((0, 1))
        dds = np.zeros((0, 3))
        sols = []

        nrestarts_requested = nrestarts
        nrestarts = 0
        for ii in range(nrestarts_requested):
            if verbose:
                oahelper.tprint('Doptim: iteration %d/%d (time %.1f/%.1f)' %
                                (ii, nrestarts, time.time() - t0, maxtime), dt=4)
            al = arrayclass.randomarray(1)

            if isinstance(optimfunc, list):
                alpha = optimfunc
                Ax = oalib.optimDeff(
                    al, arrayclass, alpha, verbose >= 2, method, niter, nabort)
                dd = Ax.Defficiencies()
                score = oalib.scoreD(dd, optimfunc)
            else:
                score, Ax = optimDeffPython(
                    al, verbose=0, niter=niter, alpha=optimfunc, method=method, nabort=nabort)
                dd = Ax.Defficiencies()

            if time.time() - t0 > maxtime:
                if verbose:
                    print('Doptim: max time exceeded, aborting')
                break

            scores = np.vstack((scores, [score]))
            dds = np.vstack((dds, dd))
            sols.append(Ax)
            nrestarts = nrestarts + 1

            if verbose >= 2:
                print('  generated array: %f %f %f' % (dd[0], dd[1], dd[2]))

            if selectpareto and ii % 502 == 0:
                scores, dds, sols = filterPareto(scores, dds, sols)

    if selectpareto:
        if verbose >= 2:
            print('Doptim: before Pareto filter (%d arrays)' % len(sols))
        scores, dds, sols = filterPareto(scores, dds, sols)

    if verbose:
        dt = time.time() - t0
        print('Doptim: done (%d arrays, %.1f [s])' % (len(sols), dt))

    # sort & select
    scores, dds, sols = selectDn(scores, dds, sols, nout=nout)

    return scores, dds, sols, nrestarts
