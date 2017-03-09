# -*- coding: utf-8 -*-
"""

Collection of helper functions for OA package

Pieter Eendebak <pieter.eendebak@gmail.com>

@author: eendebakpt
"""

#%% Load packages
from __future__ import print_function

import oapackage
import oalib
import sys
import os
import numpy as np
import functools
from collections import Counter
import operator
import fileinput
import traceback
import re
from time import gmtime, strftime
import time
import warnings


try:
    import matplotlib
    import matplotlib.pyplot as plt
except:
    warnings.warn(
        'oahelper: matplotlib cannot be found, not all functionality is available')
    pass


#%% Load Qt support

try:
    try:
        from PySide import QtGui
    except:
        from PyQt4 import QtGui

    _applocalqt = None

    def monitorSizes(verbose=0):
        """ Return monitor sizes """
        if sys.platform == 'win32':
            import ctypes
            user32 = ctypes.windll.user32
            wa = [
                [0, 0, user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)]]
        else:
            _applocalqt = QtGui.QApplication.instance()
            if _applocalqt is None:
                _applocalqt = QtGui.QApplication([])
                _qd = QtGui.QDesktopWidget()
            else:
                if 0:
                    print('get QDesktopWidget()')
                _qd = QtGui.QDesktopWidget()

            nmon = _qd.screenCount()
            wa = [_qd.screenGeometry(ii) for ii in range(nmon)]
            wa = [[w.x(), w.y(), w.width(), w.height()] for w in wa]

            if verbose:
                for ii, w in enumerate(wa):
                    print('monitor %d: %s' % (ii, str(w)))
        return wa
except:
    def monitorSizes(verbose=0):
        return [[0, 0, 1280, 720]]
    pass


def tilefigs(lst, geometry, ww=None, raisewindows=False, tofront=False, verbose=0):
    """ Tile figure windows on a specified area """
    mngr = plt.get_current_fig_manager()
    be = matplotlib.get_backend()
    if ww is None:
        ww = monitorSizes()[-1]

    w = ww[2] / geometry[0]
    h = ww[3] / geometry[1]

    # wm=plt.get_current_fig_manager()

    if verbose:
        print('tilefigs: ww %s, w %d h %d' % (str(ww), w, h))
    for ii, f in enumerate(lst):
        if not plt.fignum_exists(f):
            continue
        fig = plt.figure(f)
        iim = ii % np.prod(geometry)
        ix = iim % geometry[0]
        iy = np.floor(float(iim) / geometry[0])
        x = ww[0] + ix * w
        y = ww[1] + iy * h
        if verbose:
            print('ii %d: %d %d: f %d: %d %d %d %d' %
                  (ii, ix, iy, f, x, y, w, h))
            if verbose >= 2:
                print('  window %s' % mngr.get_window_title())
        if be == 'WXAgg':
            fig.canvas.manager.window.SetPosition((x, y))
            fig.canvas.manager.window.SetSize((w, h))
        if be == 'WX':
            fig.canvas.manager.window.SetPosition((x, y))
            fig.canvas.manager.window.SetSize((w, h))
        if be == 'agg':
            fig.canvas.manager.window.SetPosition((x, y))
            fig.canvas.manager.window.resize(w, h)
        if be == 'Qt4Agg' or be == 'QT4' or be == 'QT5Agg':
            # assume Qt canvas
            try:
                fig.canvas.manager.window.move(x, y)
                fig.canvas.manager.window.resize(w, h)
                fig.canvas.manager.window.setGeometry(x, y, w, h)
                # mngr.window.setGeometry(x,y,w,h)
            except Exception as e:
                print('problem with window manager: ', )
                print(be)
                print(e)
                pass
        if raisewindows:
            mngr.window.raise_()
        if tofront:
            plt.figure(f)


def plot2Dline(line, *args, **kwargs):
    """ Plot a 2D line in a matplotlib figure """
    if np.abs(line[1]) > .001:
        xx = plt.xlim()
        xx = np.array(xx)
        yy = (-line[2] - line[0] * xx) / line[1]
        plt.plot(xx, yy, *args, **kwargs)
    else:
        yy = np.array(plt.ylim())
        xx = (-line[2] - line[1] * yy) / line[0]
        plt.plot(xx, yy, *args, **kwargs)


#%% Make nice plots
# http://blog.olgabotvinnik.com/prettyplotlib/

def niceplot(ax, fig=None, despine=True, verbose=0, figurebg=True, tightlayout=True, legend=None, almost_black='#222222'):
    """ Create a good looking plot

    The code:
        - removes spines
        - makes legend and spines lighter
        - makes legend lighter

    """

    # Remove top and right axes lines ("spines")
    if verbose:
        print('niceplot: remove spines')

    spines_to_keep = ['bottom', 'left']
    if despine:
        spines_to_remove = ['top', 'right']
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
    else:
        spines_to_keep += ['top', 'right']
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    if verbose:
        print('niceplot: reduce spine intensity')

    for spine in spines_to_keep:
        ax.spines[spine].set_linewidth(0.5)
        ax.spines[spine].set_color(almost_black)

    ax.tick_params(axis='both', direction='out')

    if not legend is None:
        if verbose:
            print('niceplot: adjust legend')

        # Remove the line around the legend box, and instead fill it with a light grey
        # Also only use one point for the scatterplot legend because the user will
        # get the idea after just one, they don't need three.
        light_grey = np.array([float(241) / float(255)] * 3)
        rect = legend.get_frame()
        rect.set_facecolor(light_grey)
        middle_grey = np.array([float(151) / float(255)] * 3)
        rect.set_edgecolor(middle_grey)
        # rect.set_linewidth(0.0)

        # Change the legend label colors to almost black, too
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)

        # ttx=legend.get_texts()
        #[v.set_color(almost_black) for v in ttx]

    # fig.tight_layout(pad=0.1)
    if not fig is None and tightlayout:
        fig.tight_layout(pad=1.0)

    if figurebg and fig is not None:
        fig.set_facecolor('w')

    plt.draw()
    plt.show()


def enlargelims(factor=1.05):
    """ Enlarge the limits of a plot
    Example:
      >>> enlargelims(1.1)
    """
    xl = plt.xlim()
    d = (factor - 1) * (xl[1] - xl[0]) / 2
    xl = (xl[0] - d, xl[1] + d)
    plt.xlim(xl)
    yl = plt.ylim()
    d = (factor - 1) * (yl[1] - yl[0]) / 2
    yl = (yl[0] - d, yl[1] + d)
    plt.ylim(yl)

#%%


def selectParetoArrays(allarrays, pp):
    """ Select arrays using a Pareto object """
    paretoarrays = oalib.arraylist_t()
    paretoidx = np.array(pp.allindices())
    ww = oalib.longVector(tuple(paretoidx.tolist()))
    oalib.selectArrays(allarrays, ww, paretoarrays)
    return paretoarrays


def joinArrayLists(ww):
    """ Concetenate a list of array lists into a single list """
    ll = oalib.arraylist_t()
    for w in ww:
        for al in w:
            ll.push_back(al)
    return ll


def createPareto(dds, verbose=1):
    """ Create Pareto object from dataset """
    pp = oapackage.ParetoDoubleLong()

    for ii in range(dds.shape[0]):
        v = (dds[ii, :]).astype(np.float64)
        pp.addvalue(v, ii)

    if verbose:
        pp.show(verbose)
    return pp


#%% Utils


def static_var(varname, value):
    """ Helper function to create a static variable """
    def decorate(func):
        setattr(func, varname, value)
        return func
    return decorate


@static_var("time", 0)
def tprint(string, dt=1, output=False):
    """ Print progress of a loop every dt seconds """
    if (time.time() - tprint.time) > dt:
        print(string)
        tprint.time = time.time()
        if output:
            return True
        else:
            return
    else:
        if output:
            return False
        else:
            return


def timeString(tt=None):
    """ Return a string with the current time or specified time """
    if tt is None:
        tt = gmtime()
    ts = strftime("%Y-%m-%d %H:%M:%S", tt)
    return ts


def findfilesR(p, patt):
    """ Get a list of files (recursive) """
    lst = []
    for root, dirs, files in os.walk(p, topdown=False):
        lst += [os.path.join(root, f) for f in files]
        # for name in files:
        #    lst += [name]
    rr = re.compile(patt)
    lst = [l for l in lst if re.match(rr, l)]
    return lst


def findfiles(p, patt=None):
    """ Get a list of files """
    lst = os.listdir(p)
    if not patt is None:
        rr = re.compile(patt)
        lst = [l for l in lst if re.match(rr, l)]
    return lst


def finddirectories(p, patt=None):
    """ Get a list of directories """
    lst = os.listdir(p)
    if not patt is None:
        rr = re.compile(patt)
        lst = [l for l in lst if re.match(rr, l)]
    lst = [l for l in lst if os.path.isdir(os.path.join(p, l))]
    return lst


def oainfo(afile, verbose=1):
    """ Print information about a file containing arrays """
    af = oapackage.arrayfile_t(afile, verbose)
    print(af.showstr())
    af.closefile()


def oaIsBinary(afile):
    """ Return true if array file is in binary format """
    af = oapackage.arrayfile_t(afile)
    ret = af.isbinary()
    af.closefile()
    return ret


def fac(n):
    """ Return n! (factorial) """
    if n == 1 or n == 0:
        return 1
    else:
        return n * fac(n - 1)


def choose(n, k):
    """ Return n choose k """
    ntok = 1
    for t in range(min(k, n - k)):
        ntok = ntok * (n - t) // (t + 1)
    return ntok
#    return fac(n)/(fac(n-k)*fac(k))


def array2latex(X, header=1, hlines=[], floatfmt='%g', comment=None, hlinespace=None, mode='tabular', tabchar='c'):
    """ Convert numpy array to Latex tabular """
    ss = ''
    if comment is not None:
        ss += '%% %s\n' % str(comment)
    if header:
        if mode == 'tabular':
            if len(tabchar) == 1:
                cc = tabchar * X.shape[1]
            else:
                cc = tabchar + tabchar[-1] * (X.shape[1] - len(tabchar))
            ss += '\\begin{tabular}{%s}' % cc + chr(10)
        elif mode == 'psmallmatrix':
            ss += '\\begin{psmallmatrix}' + chr(10)
        else:
            ss += '\\begin{pmatrix}' + chr(10)
        # ss += '\hline' + chr(10)
    for ii in range(X.shape[0]):
        r = X[ii, :]
#        for jj in range(X.shape[1]):
#            v=r[jj]
        if isinstance(r[0], str):
            ss += ' & '.join(['%s' % x for x in r])
        else:
            ss += ' & '.join([floatfmt % x for x in r])
        if ii < (X.shape[0]) - 1 or not header:
            ss += '  \\\\' + chr(10)
        else:
            ss += '  ' + chr(10)
        if ii in hlines:
            ss += '\hline' + chr(10)
            if hlinespace != None:
                ss += '\\rule[+%.2fex]{0pt}{0pt}' % hlinespace
    if header:
        if mode == 'tabular':
            ss += '\\end{tabular}'
        elif mode == 'psmallmatrix':
            ss += '\\end{psmallmatrix}' + chr(10)
        else:
            ss += '\\end{pmatrix}' + chr(10)
    return ss


def array2latex_old(ltable, htable=None):
    """ Write a numpy array to latex format """
    print('please use array2latex from researchOA instead')
    ss = ''
    ccc = 'c' * ltable.shape[1]
    ss += '\\begin{tabular}{%s}\n' % ccc
    if not htable is None:
        ss += " \\\\\n".join(
            [" & ".join(['\\textbf{%s}' % val for val in htable])])
        ss += ' \\\\ \n'
        ss += '\\hline \n'
    for ii in range(0, ltable.shape[0]):
        line = ltable[ii, :]
        ss += " & ".join([str(val) for val in line])
        ss += ' \\\\ \n'
    ss += '\\end{tabular}\n'
    return ss


def array2html(X, header=1, tablestyle='border-collapse: collapse;', trclass='', tdstyle='', trstyle='', thstyle=''):
    """ Convert Numpy array to HTML table

    Arguments
    ---------
        X : numpy array
            array to be converted
        header : integer
            use header or not
    Returns
    -------
        page : markup html object
            generated table in HTML

    """
    page = markup.page()
    page.add('<!-- Created by array2html -->\n')
    page.table(style=tablestyle)
    offset = 0
    nc = X.shape[1]
    nr = X.shape[0]

    if isinstance(trstyle, str):
        trstyle = [trstyle] * nr
    if isinstance(trclass, str):
        trclass = [trclass] * nr

    ri = 0
    if header:
        page.tr(style='font-weight: bold; border-bottom: solid 1px black;' +
                trstyle[ri], class_=trclass[ri])
        ri = ri + 1
        for ii in range(nc):
            if isinstance(X[offset, ii], tuple):
                print('array2html: tuple instance')
                page.th(X[offset, ii][0], style=thstyle + X[offset, ii][1])
            else:
                page.th(X[offset, ii], style=thstyle)
        page.tr.close()
        offset = offset + 1

    nr = X.shape[0] - offset
    for r in range(nr):
        page.tr(style=trstyle[ri], _class=trclass[ri])
        for ii in range(nc):
            if isinstance(X[offset, ii], tuple):
                page.td(X[offset, ii][0], style=tdstyle + X[offset, ii][1])
            else:
                page.td(X[offset, ii], style=tdstyle)

        page.tr.close()
        offset = offset + 1
        ri = ri + 1
    page.table.close()
    return page

import subprocess


def runcommand(cmd, dryrun=0, idstr=None, verbose=1, logfile=None, shell=True):
    """ Run specified command in external environment 

    Returns:
        r (int): return value of the shell command
    """
    if not idstr is None:
        cmd = 'echo "idstr: %s";\n' % idstr + cmd
    if verbose >= 2:
        print('cmd: %s' % cmd)
    r = 0
    if not dryrun:

        process = subprocess.Popen(
            cmd, bufsize=1, stdout=subprocess.PIPE, shell=shell)
        for jj in range(10000000):
            r = process.poll()
            #print('poll done... %s' % r)
            line = process.stdout.readline()
            # print('poll...')
            if verbose >= 2:
                print('runcommand: jj %d' % jj)
            if verbose >= 3:
                print('runcommand: jj %d: "%s"' % (jj, line))
            if len(line) == 0:
                if verbose >= 2:
                    print('runcommand: len(line) %d' % len(line))

                break

            if r is not None:
                break
            if jj > 20000000:
                print('error: early abort of runcommand')
                break
            line = line.decode(encoding='UTF-8')
            sys.stdout.write(str(line))
            if jj % 2 == 0:
                sys.stdout.flush()
            if verbose >= 2:
                print('end of loop...')
        #print('exit loop...')
        r = process.poll()
        # r = os.system(cmd) # old method
        if (not r == 0):
            print('runcommand: cmd returned error! r=%d')
            print(cmd)
            return r
    else:
        if verbose >= 2 or (verbose and logfile is None):
            print('### dryrun\n%s\n###\n' % cmd)
    if logfile is not None:
        fid = open(logfile, 'a')
        fid.write('\n###\n')
        fid.write(cmd)
        fid.close()
    else:
        pass
        #print('no logfile')
        #raise Exception('no logfile')

    # all good
    return r


def getArrayFile(afile):
    """ Return pointer to array file
    Automatically check for compressed array files
    """
    afilegz = afile = afile + '.gz'
    if not os.path.exists(afile) and os.path.exists(afilegz):
        afile = afilegz
    else:
        afile = afile
    return afile


def checkOAfile(afile, verbose=0):
    """ Return pointer to array file
    Automatically check for compressed array files
    """
    af = oapackage.arrayfile_t(afile, verbose)
    nm = af.filename
    if af.isopen():
        return nm
    else:
        return None


def checkArrayFile(afile, cache=1):
    """ Check whether a file or list of files exists
        If the file does not exist, check whether a compressed file does exist
        cache: 0 (no), 1 (check), -1 (always)
    """
    if cache == -1:
        return True
    if cache == 0:
        return False

    if os.path.exists(afile):
        return True
    else:
        if afile.endswith('.oa'):
            if os.path.exists(afile + '.gz'):
                return True
            else:
                return False
        else:
            return False

try:
    basestring  # attempt to evaluate basestring

    def isstr(s):
        return isinstance(s, basestring)
except NameError:
    # probably Python 3.x
    def isstr(s):
        return isinstance(s, str)


def checkFiles(lst, cache=1, verbose=0):
    """ Check whether a file or list of files exists
        cache: 0 (no), 1 (check), -1 (always)

        Returns False if one or more of the files do not exist
        Returns True if all files exist
    """
    if cache == -1:
        return True
    if cache == 0:
        return False

    if isstr(lst):
        lst = [lst]
    c = True
    for f in lst:
        if not os.path.exists(f):
            if verbose:
                print('checkFiles: file %s does not exist' % f)
            c = False
            break
    return c

#%%


def checkFilesOA(lst, cache=1, verbose=0):
    """ Check whether a file or list of files exists
        cache: 0 (no), 1 (check), -1 (always)
        For array files also the .gz extension is checked

        Returns False if one or more of the files do not exist
        Returns True if all files exist
    """
    if cache == -1:
        return True
    if cache == 0:
        return False

    if isstr(lst):
        lst = [lst]
    c = True
    for f in lst:
        r = checkArrayFile(f, cache=True)
        if not r:
            if verbose:
                print('checkFiles: file %s does not exist' % f)
            c = False
            break
    return c

#%%

# import gc


def randomizearrayfile(afile, afileout, verbose=1):
    """ Randomize a file with arrays """
    lst = oapackage.readarrayfile(afile)
    rlst = oapackage.oalib.arraylist_t()
    for ii, al in enumerate(lst):
        adata = oalib.arraylink2arraydata(al)
        trans = oalib.array_transformation_t(adata)
        if verbose >= 2:
            print('randomizing array %d\n' % ii)
        trans.randomize()
        alr = trans.apply(al)
        rlst.push_back(alr)
    oapackage.writearrayfile(afileout, rlst)


def nArrayFile(afile, verbose=1):
    """ Return number of arrays in file """
    af = oalib.arrayfile_t(afile, 0)
    n = af.narrays
    af.closefile()
    return n


def selectArrays(infile, outfile, idx, afmode=oalib.ATEXT, verbose=1, cache=1):
    """ Select arrays in a file by indices

    Input:
        - infile (string)
        - outfile (string)
        - inx (list of indices)
    Output:
        - None

    """
    if not checkFiles(outfile, cache=cache):
        gidxint = oalib.intVector([int(x) for x in idx])
        # print('  xx select: %s' % str(gidxint))
        # print('  xx select: %s' % str([v for v in gidxint]))
        sols = oalib.arraylist_t()
        oalib.selectArrays(infile, gidxint, sols, 0)
        af = oalib.arrayfile_t(infile, 1)
        af.nrows
        r = oalib.writearrayfile(outfile, sols, afmode, af.nrows, af.ncols)
        if verbose >= 2:
            print('selectArrays: write array file %s in mode %d' %
                  (outfile, afmode))
    else:
        if verbose >= 2:
            print('output file %s already exists' % outfile)


def floatformat(number, mind=2, maxd=4):
    """ Format a floating point number into a string """
    rounded = int(number * (10**mind)) / float(10**mind)
    if number == rounded:
        return ('%%.0%df' % mind) % number
    else:
        return ('%%.0%dg' % maxd) % number


def safemin(data, default=0):
    """ Return mininum of array with default value for empty array """
    if data.size == 0:
        m = default
    else:
        m = data.min()
    return m


def mkdirc(xdir):
    """ Create directory """
    if not os.path.exists(xdir):
        os.mkdir(xdir)
    return xdir




def parseProcessingTime(logfile, verbose=0):
    """ Parse a log file to calculate the processing time """

    try:
        import dateutil.parser
    except:
        warnings.warn('oahelper: could not load dateutil package...')
        pass

    fileinput.close()
    tstart = None
    tend = None
    dtr = None
    try:
        for line in fileinput.input([logfile]):
            if line.startswith('#time'):
                if verbose >= 1:
                    # print(line)
                    print('parseProcessingTime: line: %s' % line, end="")
                    # print('xy %s' % line[10:])
            if line.startswith('#time start:'):
                tstart = dateutil.parser.parse(line[13:])
                if verbose >= 2:
                    print('parseProcessingTime: tstart: %s' % tstart)
            # else:
            #    print('invalid line? %s' % line)
            if line.startswith('#time end:'):
                tend = dateutil.parser.parse(line[10:])
                if verbose >= 2:
                    print('parseProcessingTime: tend: %s' % tend)
            if line.startswith('#time total:'):
                dtr = (line[13:])
                dtr = float(dtr[:-5])
                if verbose >= 2:
                    print('parseProcessingTime: total: %s' % dtr)
            # else:
            #    print('invalid line? %s' % line)
        if tstart is not None and tend is not None:
            dt = tend - tstart
            dtt = dt.total_seconds()
        else:
            dtt = -1
    except:
        if verbose:
            print('error processing log %s' % logfile)
            traceback.print_exc(file=sys.stdout)
        dtt = -1
    if not dtr is None:
        if abs(dtr - dtt) > 10:
            print(
                'parseProcessingTime: warning difference in reported times %.1f dtr %.1f [s]' % (dtt, dtr))
    return dtt


def safemax(data, default=0):
    """ Return maximum of array with default value for empty array """

    if isinstance(data, list):
        if len(data) == 0:
            m = default
        else:
            m = max(data)
        return m
    if data.size == 0:
        m = default
    else:
        m = data.max()
    return m


def series2htmlstr(ad, html=1, case=0):
    """ Convert arraydata_t to html formatted string """
#    s=[ad.s[ii] for ii in range(0, ad.ncols)]
    s = list(ad.getS())
    p = -1
    n = 0
    aa = []
    bb = []
    while n < len(s):
        if s[n] != p:
            p = s[n]
            aa.append(p)
            bb.append(1)
        else:
            bb[-1] += 1
        n = n + 1
    if case == 0:
        if bb[-1] > 1:
            bb[-1] = 'a'
    hstr = 'OA(%d; %d; ' % (ad.N, ad.strength)
    for ii in range(0, len(aa)):
        if html:
            hstr += '%d<sup>%s</sup>' % (aa[ii], str(bb[ii]))
        else:
            hstr += '%d^%s' % (aa[ii], str(bb[ii]))
    hstr += ')'
    return hstr


def gwlp2str(gmadata, t=None, sformat=None, jstr=','):
    """ Convert GWLP value to string format """
    if gmadata is None:
        return '-'
    if isinstance(gmadata, tuple):
        # do nothing
        gmadata
    else:
        if isinstance(gmadata, list):
            # do nothing
            gmadata
        else:
            gmadata[gmadata < 0] = 0
        if not(np.abs(gmadata[0] - 1) < 1e-12 and np.abs(gmadata[1]) < 1e-12):
            print('warning: data are not good GWPL data!!!!')
            return ''
            # pdb.set_trace()
    bgma = np.around(gmadata, decimals=12)
    if not t is None:
        bgma = bgma[(t + 1):]
    if sformat is None:
        gstr = jstr.join([floatformat(v, mind=2, maxd=4) for v in bgma])
    else:
        gstr = jstr.join([sformat % v for v in bgma])
    return gstr


def selectJ(sols0, jj=5, jresults=None, verbose=1):
    """ Select only arrays with J-characteristics non-zero """
    if jresults is None:
        jresults = oalib.analyseArrays(sols0, verbose, jj)

    solseo = oalib.arraylist_t()
    v = []
    for ii, js in enumerate(jresults):
        # js.show()
        #    if not np.array(js.vals)[0]==0:
        #        solseo.append(sols0[ii])
        v.append(js.vals[0])

    si = [i for (i, j) in sorted(enumerate(v), key=operator.itemgetter(1))]
    for jj in si:
        if v[jj] > 0:
            solseo.append(sols0[jj])
    if verbose:
        print('selectJ: kept %d/%d solutions' % (solseo.size(), sols0.size()))
    return solseo


def extendSingleArray(A, adata, t=3, verbose=1):
    """ Extend a single array """
    oaoptions = oalib.OAextend()
    sols0 = oalib.arraylist_t()
    solsx = oalib.arraylist_t()
    sols0.push_back(A)
    #N = A.n_rows
    k = A.n_columns
    n = oalib.extend_arraylist(sols0, adata, oaoptions, k, solsx)
    sys.stdout.flush()
    return solsx


def runExtend(N, k, t=3, l=2, verbose=1, initsols=None, nums=[], algorithm=None):
    """ Run extension algorithm and return arrays

    Arguments
    ---------
    N : integer
        number of rows
    k: integer
        number of columns
    t: integer
        strength of the arrays

    >>> r = runExtend(16, 5, 3, verbose=0)    
    """
    if verbose:
        print('runExtend: N=%d, k=%d, t=%d' % (N, k, t))
    if isinstance(l, list):  # types.ListType):
        ll = l
    else:
        ll = [l]
    ll = ll + [ll[-1]] * (k - len(ll))
    s = oalib.intVector(ll)
    adata = oalib.arraydata_t(s, N, t, k)
    al = oalib.array_link(adata.N, adata.strength, 1)
    al.create_root(adata)
    if initsols is None:
        sols0 = oalib.arraylist_t()
        sols0.append(al)
        tstart = t
    else:
        sols0 = initsols
        tstart = sols0[0].n_columns

    oaoptions = oalib.OAextend()
    if algorithm is None:
        oaoptions.setAlgorithmAuto(adata)
    else:
        oaoptions.setAlgorithm(algorithm, adata)
    solsx = sols0
    for ii in range(tstart, k):
        solsx = oalib.arraylist_t()
        oalib.extend_arraylist(sols0, adata, oaoptions, ii, solsx)
        if verbose >= 2:
            print(' ii %d: %d' % (ii, solsx.size()))
        sols0 = solsx
        nums.append(solsx.size())
        sys.stdout.flush()
    return solsx


def compressOAfile(afile, decompress=False, verbose=1):
    """ Compress an OA array file """
    af = oalib.arrayfile_t(afile, 0)
    if verbose >= 2:
        print('file %s: binary %s' % (afile, af.isbinary()))
    if not sys.platform == 'linux2':
        if verbose:
            print('compressOAfile: not compressing file %s (platform not supported)' %
                  afile)
        return False
    if af.isbinary() and not af.iscompressed and sys.platform == 'linux2':
        if verbose:
            print('compressOAfile: compressing file %s' % afile)
        cmd = 'gzip -f %s' % afile
        r = os.system(cmd)
        if not r == 0:
            print('compressOAfile: error compressing file %s' % afile)
        return True
    else:
        if not af.isbinary():
            if verbose:
                print('compressOAfile: not compressing file %s (file is in text mode)' %
                      afile)
        else:
            if verbose:
                print('compressOAfile: not compressing file %s (file %s is compressed)' %
                      (afile, af.filename))
        return False

if 0:
    def pointer2np(p, sz):
        """ Convert array pointer to numpy array """
        if isinstance(sz, tuple):
            1
        else:
            sz = (sz,)
        ip = oalib.intArray.frompointer(p)
        A = np.zeros(sz)
        for ii in range(0, np.prod(sz)):
            A[ii] = ip[ii]
        A = A.reshape(sz)
        return A

if 0:
    import ctypes

    def arraylink2array(al):
        """ Convert arraylink to numpy array """
        print('use araylink.getarray')
        w = al.array
        nrows = al.n_rows
        ncols = al.n_columns
        func = ctypes.pythonapi.PyBuffer_FromMemory
        func.restype = ctypes.py_object
        buffer = func(w.__long__(), nrows * ncols * 2 * 2)
        xx = np.frombuffer(buffer, np.int32, count=nrows * ncols)
        xx = xx.reshape((ncols, nrows)).transpose()
        return xx


def getposjvals(A, t, verbose=1):
    N = A.shape[0]
    k = A.shape[1]
    ncols = t + 1
    jstep = 2**(t + 1)
    njvals = int(1 + (float(N) / jstep))
    jvals = [N - jstep * l for l in range(0, njvals)]
    return jvals


def arraystats(A, verbose=1):
    """ Return statistics of an array """
    Af = A.transpose().flatten()
    al = oalib.array_link(A.shape[0], A.shape[1], 0)
    alist = oalib.arraylist_t(1)
    alist[0] = al
    # al.array
    ia = oalib.intArray.frompointer(al.array)
    for x in range(0, A.size):
        ia[x] = int(Af[x])
#    al.showarray()
    jresults = oalib.analyseArrays(alist, 0)
    js = jresults[0]
    vals = pointer2np(js.vals, js.nc)
    jvals = getposjvals(A, 3, verbose=0)

    jv = np.zeros(len(jvals))
    c = Counter(abs(vals))
    for ii, xx in enumerate(jvals):
        if c.has_key(xx):
            jv[ii] = c.get(xx)
    if verbose:
        print('Possible j-values: %s' % jvals, end="")
        print('     values: %s' % jv.astype(int))

    N = A.shape[0]
    Ak = (1 / N ^ 2) * sum(vals**2)
    print('Ak: %s' % Ak)
    return Ak


def argsort(seq):
    """ Stable argsort """
    # http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    return sorted(range(len(seq)), key=seq.__getitem__)


def jseq(xx, comb):
    pp = functools.reduce(lambda y, i: xx[:, i] + y, comb, 0)
    jseq = 1 - 2 * np.mod(pp, 2)
    return jseq


def sortrows(x):
    """ Sort rows of an array, return indices"""
    if len(x.shape) == 1:
        nn = 1
        sind = np.argsort(x)  # or np.argsort(x, order=('x', 'y'))
        return sind
    else:
        nn = x.shape[1]
    dt = [('x%d' % xx, 'f8') for xx in range(0, nn)]
    if x.size == 0:
        # hack
        sind = np.array([])
        return sind
    v = [tuple(x[kk, :]) for kk in range(0, x.shape[0])]
    w = np.array(v, dtype=dt)
    sind = np.argsort(w)  # or np.argsort(x, order=('x', 'y'))
    return sind


def sortcols(X):
    """ Sort columns of an array, return indices"""
    sind = sortrows(X.transpose())
    return sind


def showtriangles(jresults, showindex=1):
    """ Show triangle of j-values """
    if isinstance(jresults, oalib.jstruct_t):
        showindex = 0
        jresults = (jresults,)

    js = jresults[0]
    i = 0
    idx = [i]
    for j in range(4, js.k + 1):
        nn = choose(j - 1, 3)
        i = i + nn
        idx.append(i)
    for jj, js in enumerate(jresults):
        vals = oalib.intArray.frompointer(js.vals)
        if showindex:
            print('i: %d' % jj)
        for kk in range(0, js.k - 3):
            xx = [vals[v] for v in range(idx[kk], idx[kk + 1])]
            s = ','.join(map(str, xx))
            print('%s' % s)
#%%

from oapackage import markup
import webbrowser
import tempfile


def testHtml(hh=None):
    """ Test a short snippet of HTML """
    if hh is None:
        return
    page = markup.page()
    page.init()
    page.body()
    page.add(hh)
    page.body.close()
    tmp, f = tempfile.mkstemp('.html')
    with open(f, 'wt') as fname:
        fname.write(str(page))
        fname.close()
        webbrowser.open(fname.name)


def designStandardError(al):
    """ Return standard errors for a design

    Arguments
    ---------
    al : array
            design

    Output
    ------
    m0, m1, m2 : arrays
            standard errors

    """

    X = np.matrix(al.getModelMatrix(2))
    k = al.n_columns

    m = 1 + k + k * (k - 1) / 2

    # scalefac=np.sqrt(al.n_rows)
    scalefac = 1
    M = (X.transpose() * X / scalefac).I

    mm = np.array(M.diagonal()).flatten()

    m1 = mm[1:(1 + k)].flatten()
    m1 = m1[np.argsort(m1)]
    m2 = mm[(1 + k):]
    m2 = m2[np.argsort(m2)]
    m0 = mm[0]
    return np.sqrt(m0), np.sqrt(m1), np.sqrt(m2)
    # return m0,m1, m2


#%%
def DefficiencyBound(D, k, k2):
    """ Calculate the D-efficiency bound of an array extension

    Input
    -----
    D : float
        D-efficiency of the design
    k, k2: integers
        numbers of columns 
    Output
    ------
    D2 : float
        bound on the D-efficiency of extensions of a design with k columns to k2 columns    

    """
    m = 1. + k + k * (k - 1) / 2
    m2 = 1. + k2 + k2 * (k2 - 1) / 2
    D2 = D**(m / m2)
    return D2

#%% Misc


def setWindowRectangle(x, y=None, w=None, h=None, mngr=None, be=None):
    """ Position the current Matplotlib figure at the specified position
    Usage: setWindowRectangle(x,y,w,h)
    """
    if y is None:
        y = x[1]
        w = x[2]
        h = x[3]
        x = x[0]
    if mngr is None:
        mngr = plt.get_current_fig_manager()
    be = matplotlib.get_backend()
    if be == 'WXAgg':
        mngr.canvas.manager.window.SetPosition((x, y))
        mngr.canvas.manager.window.SetSize((w, h))
    elif be == 'agg':
            # mngr.canvas.manager.window.setGeometry(x,y,w,h)
        mngr.canvas.manager.window.SetPosition((x, y))
        mngr.canvas.manager.window.resize(w, h)
    elif be == 'module://IPython.kernel.zmq.pylab.backend_inline':
        pass
        # mngr.canvas.manager.window.setGeometry(x,y,w,h)
        #mngr.canvas.manager.window.SetPosition((x, y))
        #mngr.canvas.manager.window.resize(w, h)
    else:
        # assume Qt canvas
        mngr.canvas.manager.window.move(x, y)
        mngr.canvas.manager.window.resize(w, h)
        mngr.canvas.manager.window.setGeometry(x, y, w, h)
        # mngr.window.setGeometry(x,y,w,h)


def makearraylink(al):
    """ Convert array to array_link object """
    if isinstance(al, np.ndarray):
        tmp = oalib.array_link()
        tmp.setarray(al)
        al = tmp
    return al


def formatC(al, wrap=True):
    """ Format array the inclusion in C code """
    l = np.array(al).T.flatten().tolist()
    s = ','.join(['%d' % x for x in l])
    if wrap:
        s = '\tarray_link al ( %d,%d, 0 );\n\tint tmp[] = {' % (al.n_rows, al.n_columns) + s + '};'
    return s

#%%


def create_pareto_element(values, pareto=None):
    """ Create a vector of mvalue_t elements 
    Args:
        vv (list): list with tuples or arrays
    """
    if pareto is None:
        pareto = oalib.ParetoDoubleLong()
    if isinstance(pareto, oalib.ParetoMultiLongLong):
        vector_pareto = oalib.mvalueVector()
        for v in values:
            vec = oalib.mvalue_t_long(list(v))
            vector_pareto.push_back(vec)
    else:
        vector_pareto = vv
    return vector_pareto
