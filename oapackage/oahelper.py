# -*- coding: utf-8 -*-
""" Collection of helper functions for OA package

@author: Pieter Eendebak <pieter.eendebak@gmail.com>
"""

# %% Load packages
from __future__ import print_function

import sys
import os
import logging
import numpy as np
import functools
import operator
import inspect
import fileinput
import traceback
import re
from time import gmtime, strftime
import time
import warnings
import webbrowser
import tempfile
import subprocess
import dateutil.parser

try:
    import matplotlib
    import matplotlib.pyplot as plt
except BaseException:
    warnings.warn(
        'oahelper: matplotlib cannot be found, not all functionality is available')

import oapackage
import oalib


from oapackage import markup


def deprecated(func):
    """ This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used. """

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        try:
            filename = inspect.getfile(func)
        except BaseException:
            filename = '?'
        try:
            lineno = inspect.getlineno(func)
        except BaseException:
            lineno = -1
        warnings.warn_explicit(
            "Call to deprecated function {}.".format(func.__name__),
            category=UserWarning,
            filename=filename,
            lineno=lineno,
        )
        return func(*args, **kwargs)
    return new_func

# %%


try:
    try:
        import qtpy.QtGui
        import qtpy.QtWidgets
    except BaseException:
        # no Qt support
        pass

    _applocalqt = None

    def monitorSizes(verbose=0):
        """ Return monitor sizes """
        if sys.platform == 'win32':
            import ctypes
            user32 = ctypes.windll.user32
            wa = [
                [0, 0, user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)]]
        else:
            _applocalqt = qtpy.QtWidgets.QApplication.instance()
            if _applocalqt is None:
                _applocalqt = qtpy.QtWidgets.QApplication([])
                _qd = qtpy.QtWidgets.QDesktopWidget()
            else:
                _qd = qtpy.QtWidgets.QDesktopWidget()

            nmon = _qd.screenCount()
            wa = [_qd.screenGeometry(ii) for ii in range(nmon)]
            wa = [[w.x(), w.y(), w.width(), w.height()] for w in wa]

            if verbose:
                for ii, w in enumerate(wa):
                    print('monitor %d: %s' % (ii, str(w)))
        return wa
except BaseException:
    def monitorSizes(verbose=0):
        """ Dummy implementation """
        return [[0, 0, 1280, 720]]


def guitest_monitorSizes():
    sizes = monitorSizes()
    assert(isinstance(sizes, list))


def tilefigs(lst, geometry, ww=None, raisewindows=False, tofront=False, verbose=0):
    """ Tile figure windows on a specified area """
    mngr = plt.get_current_fig_manager()
    be = matplotlib.get_backend()
    if ww is None:
        ww = monitorSizes()[-1]

    w = ww[2] / geometry[0]
    h = ww[3] / geometry[1]

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
            except Exception as ex:
                print('problem with window manager: ', )
                print('backend %s' % (be,))
                logging.exception(ex)
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


# %% Make nice plots

def niceplot(ax, fig=None, despine=True, verbose=0, figurebg=True,
             tightlayout=True, legend=None, almost_black='#222222'):
    """ Create a good looking plot

    The code:
        - removes spines
        - makes legend and spines lighter
        - makes legend lighter

    """

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

    if legend is not None:
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

    if (fig is not None) and tightlayout:
        fig.tight_layout(pad=1.0)

    if figurebg and fig is not None:
        fig.set_facecolor('w')

    plt.draw()
    plt.show()


def enlargelims(factor=1.05):
    """ Enlarge the limits of a plot

    Args:
        factor (float): factor by which to make the plot margins wider
    """
    xl = plt.xlim()
    d = (factor - 1) * (xl[1] - xl[0]) / 2
    xl = (xl[0] - d, xl[1] + d)
    plt.xlim(xl)
    yl = plt.ylim()
    d = (factor - 1) * (yl[1] - yl[0]) / 2
    yl = (yl[0] - d, yl[1] + d)
    plt.ylim(yl)

# %%


def helmert_contrasts(number_of_levels, verbose=0):
    """ Calculate Helmert contrasts for a given number of levels

    The Helmert contrasts are orthogonal and normalize such that the square equals to number of levels.
    Args:
        number_of_levels: number of levels in the number_of_levels
    Returns:
        array: array with calculated Helmert contrasts
    """
    N = number_of_levels
    meoffset = 0
    md = number_of_levels - 1

    main_effects = np.zeros((number_of_levels, md))
    Z = np.zeros((number_of_levels, md + 1))

    for value in range(number_of_levels):
        for ii in range(0, md + 1):
            Z[value, ii] = value > ii - 1

        # make Helmert contrasts (these are automatically orthogonal)
        Z[value, 0] = 1
        if (value > 0):
            Z[value, value] = value

        for q in range(1, value):
            Z[value, q] = 0
        for q in range(value + 1, md + 1):
            Z[value, q] = -1

    if verbose:
        print('helmert_contrasts: %d\n' % (number_of_levels))
        print('Z (initial creation)')
        print(Z)

    # normalize the contrasts
    for ii in range(0, md):
        normalization = Z[:, (ii + 1)].T.dot(Z[:, ii + 1])
        if verbose:
            print('helmert_contrasts: normalize number_of_levels tmp: %s ' % (normalization,))
        main_effects[:, meoffset + ii:(meoffset + ii + 1)] = np.sqrt((N)) * \
            Z[:, (ii + 1):(ii + 2)] / np.sqrt((normalization))

    return main_effects


def selectParetoArrays(array_list, pareto_object):
    """ Select arrays using a Pareto object

    Args:
        array_list (list): list of arrays
        pareto (object): oapackage Pareto object
    Returns:
        list: list with all Pareto optimal designs
    """
    paretoarrays = oalib.arraylist_t()
    paretoidx = np.array(pareto_object.allindices())
    ww = oalib.longVector(tuple(paretoidx.tolist()))
    oalib.selectArrays(array_list, ww, paretoarrays)
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


# %% Utils


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
    """ Return a string with the current time or specified time

    Args:
        tt (struct_time or None): time to convert
    Returns:
        str: formatted time
    """
    if tt is None:
        tt = gmtime()
    ts = strftime("%Y-%m-%d %H:%M:%S", tt)
    return ts


def findfilesR(directory, pattern):
    """ Get a list of files (recursive)

    Args:
        directory (str): directory to search
        patt (str): pattern
    Returns:
        list: list of matched files
    """
    lst = []
    for root, _, files in os.walk(directory, topdown=False):
        lst += [os.path.join(root, f) for f in files]
    rr = re.compile(pattern)
    lst = [l for l in lst if re.match(rr, l)]
    return lst


def findfiles(directory, pattern=None):
    """ Get a list of files
    Args:
        directory (str): directory to search
        patt (str): pattern
    Returns:
        list: list of matched files
    """
    lst = os.listdir(directory)
    if pattern is not None:
        rr = re.compile(pattern)
        lst = [l for l in lst if re.match(rr, l)]
    return lst


def finddirectories(directory, pattern=None):
    """ Get a list of directories
    Args:
        directory (str): directory to search
        patt (str): pattern
    Returns:
        list: list of matched directories
    """
    lst = os.listdir(directory)
    if pattern is not None:
        rr = re.compile(pattern)
        lst = [l for l in lst if re.match(rr, l)]
    lst = [l for l in lst if os.path.isdir(os.path.join(directory, l))]
    return lst


def oainfo(filename, verbose=1):
    """ Print information about a file containing arrays """
    af = oapackage.arrayfile_t(filename, verbose)
    print(af.showstr())
    af.closefile()


def oaIsBinary(filename):
    """ Return true if array file is in binary format """
    af = oapackage.arrayfile_t(filename)
    ret = af.isbinary()
    af.closefile()
    return ret


def fac(n):
    """ Return n! (factorial)

    Args:
        n (int): number to calculate factorial
    Returns
        integer: factorial of number
    """
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


def array2latex(X, header=1, hlines=(), floatfmt='%g', comment=None, hlinespace=None, mode='tabular', tabchar='c'):
    """ Convert numpy array to Latex tabular """
    ss = ''
    if comment is not None:
        if isinstance(comment, list):
            for line in comment:
                ss += '%% %s\n' % str(line)
        else:
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
    for ii in range(X.shape[0]):
        r = X[ii, :]
        if isinstance(r[0], str):
            ss += ' & '.join(['%s' % x for x in r])
        else:
            ss += ' & '.join([floatfmt % x for x in r])
        if ii < (X.shape[0]) - 1 or not header:
            ss += '  \\\\' + chr(10)
        else:
            ss += '  ' + chr(10)
        if ii in hlines:
            ss += r'\hline' + chr(10)
            if hlinespace is not None:
                ss += '\\rule[+%.2fex]{0pt}{0pt}' % hlinespace
    if header:
        if mode == 'tabular':
            ss += '\\end{tabular}'
        elif mode == 'psmallmatrix':
            ss += '\\end{psmallmatrix}' + chr(10)
        else:
            ss += '\\end{pmatrix}' + chr(10)
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
    for _ in range(nr):
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


def runcommand(cmd, dryrun=0, idstr=None, verbose=1, logfile=None, shell=True):
    """ Run specified command in external environment

    Returns:
        r (int): return value of the shell command
    """
    if idstr is not None:
        cmd = 'echo "idstr: %s";\n' % idstr + cmd
    if verbose >= 2:
        print('cmd: %s' % cmd)
    r = 0
    if not dryrun:

        process = subprocess.Popen(
            cmd, bufsize=1, stdout=subprocess.PIPE, shell=shell)
        for jj in range(10000000):
            r = process.poll()
            line = process.stdout.readline()
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
        r = process.poll()
        if (not r == 0):
            print('runcommand: cmd returned error! r=%d' % str(r))
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

    Args:
        afile (str): location of array file
    Returns:
        arrayfile_t or None
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
    # attempt to evaluate basestring
    basestring   # type: ignore

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


def test_checkFiles():
    def touch(fname):
        if os.path.exists(fname):
            os.utime(fname, None)
        else:
            open(fname, 'a').close()

    lst = [tempfile.mktemp()]
    r = checkFiles(lst, cache=1, verbose=1)
    assert(r is False)
    touch(lst[0])
    r = checkFiles(lst, cache=1, verbose=1)
    assert(r is True)

# %%


def checkFilesOA(lst, cache=1, verbose=0):
    """ Check whether a file or list of files exists

        Args:
            lst (list): list of files
            cache (int): 0 (always return False), 1 (check), -1 (always return True)

        For array files also the .gz extension is checked
        Returns False if one or more of the files do not exist
        Returns True if all files exist
    """
    if verbose >= 2:
        print('checkFilesOA: cache %s' % (cache, ))
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

# %%


def randomizearrayfile(input_filename, output_filename, verbose=1):
    """ Randomize a file with orthogonal arrays

    Each array is transformed with a random transformation

    Args:
        afile (str): input file
        afileout (str): output file
    """
    lst = oapackage.readarrayfile(input_filename)
    rlst = oapackage.oalib.arraylist_t()
    for ii, al in enumerate(lst):
        adata = oalib.arraylink2arraydata(al)
        trans = oalib.array_transformation_t(adata)
        if verbose >= 2:
            print('randomizing array %d\n' % ii)
        trans.randomize()
        alr = trans.apply(al)
        rlst.push_back(alr)
    oapackage.writearrayfile(output_filename, rlst)


def nArrayFile(filename, verbose=1):
    """ Return number of arrays in file

    Args:
        afile (str): name of array file
    """
    af = oalib.arrayfile_t(filename, 0)
    n = af.narrays
    af.closefile()
    return n


def selectArraysInFile(infile, outfile, idx, afmode=oalib.ATEXT, verbose=1, cache=1):
    """ Select arrays in a file by indices

    Args:
        infile (str): file with designs
        outfile (str): output  file with designs
        inx (list of indices)
    Returns:
        None

    """
    if not checkFiles(outfile, cache=cache):
        gidxint = oalib.intVector([int(x) for x in idx])
        sols = oalib.arraylist_t()
        oalib.selectArrays(infile, gidxint, sols, 0)
        af = oalib.arrayfile_t(infile, 1)
        _ = oalib.writearrayfile(outfile, sols, afmode, af.nrows, af.ncols)
        if verbose >= 2:
            print('selectArrays: write array file %s in mode %d' %
                  (outfile, afmode))
    else:
        if verbose >= 2:
            print('output file %s already exists' % outfile)

selectArrays = deprecated(selectArraysInFile)


def floatformat(number, mind=2, maxd=4):
    """ Format a floating point number into a string """
    rounded = int(number * (10**mind)) / float(10**mind)
    if number == rounded:
        return ('%%.0%df' % mind) % number
    else:
        return ('%%.0%dg' % maxd) % number


def safemin(data, default=0):
    """ Return mininum of array with default value for empty array

    Args:
        data (array or list): data to return the maximum
        default (obj): default value
    Returns:
        object: minimum value in the array or the default value

    """
    if isinstance(data, list):
        if len(data) == 0:
            minimum_value = default
        else:
            minimum_value = min(data)
        return minimum_value
    if data.size == 0:
        minimum_value = default
    else:
        minimum_value = data.min()
    return minimum_value


def safemax(data, default=0):
    """ Return maximum of array with default value for empty array

    Args:
        data (array or list ): data to return the maximum
        default (obj): default value
    Returns:
        object: maximum value in the array or the default value
    """

    if isinstance(data, list):
        if len(data) == 0:
            maximum_value = default
        else:
            maximum_value = max(data)
        return maximum_value
    if data.size == 0:
        maximum_value = default
    else:
        maximum_value = data.max()
    return maximum_value


def mkdirc(directory_name):
    """ Create directory """
    if not os.path.exists(directory_name):
        os.mkdir(directory_name)
    return directory_name


def parseProcessingTime(logfile, verbose=0):
    """ Parse a log file to calculate the processing time """

    import dateutil.parser

    fileinput.close()
    tstart = None
    tend = None
    dtr = None
    try:
        for line in fileinput.input([logfile]):
            if line.startswith('#time'):
                if verbose >= 1:
                    print('parseProcessingTime: line: %s' % line, end="")
            if line.startswith('#time start:'):
                tstart = dateutil.parser.parse(line[13:])
                if verbose >= 2:
                    print('parseProcessingTime: tstart: %s' % tstart)
            elif line.startswith('#time end:'):
                tend = dateutil.parser.parse(line[10:])
                if verbose >= 2:
                    print('parseProcessingTime: tend: %s' % tend)
            elif line.startswith('#time total:'):
                dtr = (line[13:])
                dtr = float(dtr[:-5])
                if verbose >= 2:
                    print('parseProcessingTime: total: %s' % dtr)
            else:
                pass
        if tstart is not None and tend is not None:
            dt = tend - tstart
            dtt = dt.total_seconds()
        else:
            dtt = -1
    except BaseException:
        if verbose:
            print('error processing log %s' % logfile)
            traceback.print_exc(file=sys.stdout)
        dtt = -1
    if dtr is not None:
        if abs(dtr - dtt) > 10:
            print(
                'parseProcessingTime: warning difference in reported times %.1f dtr %.1f [s]' % (dtt, dtr))
    return dtt


def series2htmlstr(ad, html=1, case=0):
    """ Convert arraydata_t to html formatted string """
    s = list(ad.factor_levels())
    p = -1
    n = 0
    levels = []
    bb = []
    while n < len(s):
        if s[n] != p:
            p = s[n]
            levels.append(p)
            bb.append(1)
        else:
            bb[-1] += 1
        n = n + 1
    if case == 0:
        if bb[-1] > 1:
            bb[-1] = 'a'
    hstr = 'OA(%d; %d; ' % (ad.N, ad.strength)
    for ii, _ in enumerate(levels):
        if html:
            hstr += '%d<sup>%s</sup>' % (levels[ii], str(bb[ii]))
        else:
            hstr += '%d^%s' % (levels[ii], str(bb[ii]))
    hstr += ')'
    return hstr


def gwlp2str(gmadata, t=None, sformat=None, jstr=','):
    """ Convert GWLP value to string format """
    if gmadata is None:
        return '-'
    if isinstance(gmadata, tuple):
        pass
    else:
        if isinstance(gmadata, list):
            pass
        else:
            gmadata[gmadata < 0] = 0
        if not(np.abs(gmadata[0] - 1) < 1e-12 and np.abs(gmadata[1]) < 1e-12):
            warnings.warn('data does not represent GWPL data', UserWarning)
            return ''
    bgma = np.around(gmadata, decimals=12)
    if t is not None:
        bgma = bgma[(t + 1):]
    if sformat is None:
        gstr = jstr.join([floatformat(v, mind=2, maxd=4) for v in bgma])
    else:
        gstr = jstr.join([sformat % v for v in bgma])
    return gstr


def selectJ(sols0, jj=5, jresults=None, verbose=1):
    """ Select only arrays with J-characteristics non-zero

    We asssume the designs are in even-odd ordering (i.e. only check the J value of first columns)
    """
    if jresults is None:
        jresults = oalib.analyseArrays(sols0, verbose, jj)

    solseo = oalib.arraylist_t()
    v = []
    for ii, js in enumerate(jresults):
        v.append(js.values[0])

    si = [i for (i, j) in sorted(enumerate(v), key=operator.itemgetter(1))]
    for jj in si:
        if v[jj] > 0:
            solseo.append(sols0[jj])
    if verbose:
        print('selectJ: kept %d/%d solutions' % (solseo.size(), len(sols0)))
    return solseo


def test_selectJ():
    al = oapackage.exampleArray(12, 1)
    sols0 = [al]
    r = selectJ(sols0)
    assert(len(r) == 0)


def extendSingleArray(A, adata, t=3, verbose=1):
    """ Extend a single orthogonal array """
    oaoptions = oalib.OAextend()
    oaoptions.setAlgorithmAuto(adata)
    sols0 = oalib.arraylist_t()
    solsx = oalib.arraylist_t()
    sols0.push_back(A)
    k = A.n_columns
    n = oalib.extend_arraylist(sols0, adata, oaoptions, k, solsx)
    assert(n >= len(solsx))
    sys.stdout.flush()
    return solsx


def test_extendSingleArray():
    A = oapackage.exampleArray(4, 1)
    adata = oapackage.arraylink2arraydata(A, extracols=2)
    B = A.selectFirstColumns(5)
    ee = extendSingleArray(B, adata, t=2, verbose=1)
    assert(ee[0].n_columns == B.n_columns + 1)
    assert(ee[1] == A.selectFirstColumns(6))


def runExtend(N, k, t=3, l=2, verbose=1, initsols=None, nums=[], algorithm=None):
    """ Run extension algorithm and return arrays

    Args:
      N (int): number of rows
      k (int): number of columns to extend to
      t (int): strength of the arrays
      l (int): factors of the designs
      initsols (None or list): list of arrays to extend, None to start with root

    Returns:
        list: list of generated designs

    Example:
       >>> import oapackage
       >>> designs = oapackage.oahelper.runExtend(16, 5, 3, verbose=0)
    """
    if verbose:
        print('runExtend: N=%d, k=%d, t=%d' % (N, k, t))
    if isinstance(l, list):
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
    """ Compress an OA array file

    Args:
        afile (str): array to compress
        decompress (bool): If True then decompress
        verbose (int): verbosity level
    """
    af = oalib.arrayfile_t(afile, 0)
    if decompress:
        raise NotImplementedError('decompressing file not implemted')

    if verbose >= 2:
        print('file %s: binary %s' % (afile, af.isbinary()))
    if not (sys.platform == 'linux2' or sys.platform == 'linux'):
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


def argsort(seq):
    """ Stable version of argsort """
    # http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    return sorted(range(len(seq)), key=seq.__getitem__)


@deprecated
def jseq(xx, comb):
    pp = functools.reduce(lambda y, i: xx[:, i] + y, comb, 0)
    jseq = 1 - 2 * np.mod(pp, 2)
    return jseq


def sortrows(x):
    """ Sort rows of an array, return indices

    The array is sorted low to high with the first column as the major column.

    Args:
        array (numpy array)
    Returns:
        list: indices of rows of sorted array
    """
    if len(x.shape) == 1:
        nn = 1
        sind = np.argsort(x)
        return sind
    else:
        nn = x.shape[1]
    dt = [('x%d' % xx, 'f8') for xx in range(0, nn)]
    if x.size == 0:
        sind = np.array([])
        return sind
    v = [tuple(x[kk, :]) for kk in range(0, x.shape[0])]
    w = np.array(v, dtype=dt)
    sind = np.argsort(w)
    return sind


def sortcols(array):
    """ Sort columns of an array, return indices

    The array is sorted low to high with the row column as the major row.

    Args:
        array (numpy array)
    Returns:
        list: indices of columns of sorted array
    """
    sind = sortrows(array.transpose())
    return sind

# %%


def testHtml(html_code=None):
    """ Test a short snippet of HTML """
    if html_code is None:
        return
    page = markup.page()
    page.init()
    page.body()
    page.add(html_code)
    page.body.close()
    _, f = tempfile.mkstemp('.html')
    with open(f, 'wt') as fname:
        fname.write(str(page))
        fname.close()
        webbrowser.open(fname.name)


def designStandardError(al):
    """ Return standard errors for a design

    Args:
      al (array): design

    Returns:
        array: array with standard errors

    """

    X = np.array(al.getModelMatrix(2))
    k = al.n_columns

    scalefac = 1
    M = np.linalg.inv(X.transpose().dot(X) / scalefac)

    mm = np.array(M.diagonal()).flatten()

    m1 = mm[1:(1 + k)].flatten()
    m1 = m1[np.argsort(m1)]
    m2 = mm[(1 + k):]
    m2 = m2[np.argsort(m2)]
    m0 = mm[0]
    return np.sqrt(m0), np.sqrt(m1), np.sqrt(m2)


# %%
def DefficiencyBound(D, k, k2):
    """ Calculate the D-efficiency bound of an array extension

    Args:
        D (float): D-efficiency of the design
        k (int): numbers of columns
        k2 (int): numbers of columns

    Returns:
        float: bound on the D-efficiency of extensions of a design with k columns to k2 columns

    """
    m = 1. + k + k * (k - 1) / 2
    m2 = 1. + k2 + k2 * (k2 - 1) / 2
    Dbound = D**(m / m2)
    return Dbound

# %% Misc


def setWindowRectangle(x, y=None, w=None, h=None, mngr=None, be=None):
    """ Position the current Matplotlib figure at the specified position

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
        mngr.canvas.manager.window.SetPosition((x, y))
        mngr.canvas.manager.window.resize(w, h)
    elif be == 'module://IPython.kernel.zmq.pylab.backend_inline':
        pass
    else:
        # assume Qt canvas
        mngr.canvas.manager.window.move(x, y)
        mngr.canvas.manager.window.resize(w, h)
        mngr.canvas.manager.window.setGeometry(x, y, w, h)


def makearraylink(array):
    """ Convert array to array_link object

    Args:
        array (numpy array): array to convert
    Returns:
        array_link
    """
    if isinstance(array, np.ndarray):
        tmp = oalib.array_link()
        tmp.setarray(array)
        array = tmp
    return array


def formatC(al, wrap=True):
    """ Format array for inclusion in C code """
    l = np.array(al).T.flatten().tolist()
    s = ','.join(['%d' % x for x in l])
    if wrap:
        s = '\tarray_link array ( %d,%d, 0 );\n\tint array_data_tmp[] = {' % (
            al.n_rows, al.n_columns) + s + '};'
        s += '\tarray.setarraydata(array_data_tmp, array.n_rows * array.n_columns);\n'
    return s


def create_pareto_element(values, pareto=None):
    """ Create a vector of mvalue_t elements
    Args:
        vv (list): list with tuples or arrays
    Returns:
        object: Pareto element
    """
    if pareto is None:
        pareto = oalib.ParetoDoubleLong()
    if isinstance(pareto, oalib.ParetoMultiLongLong):
        vector_pareto = oalib.mvalueVector()
        for v in values:
            if isinstance(v, (int)):
                # convert to list type
                v = [v]
            if not isinstance(v, (list, type)):
                raise Exception('creating Pareto element for Pareto object of type %s and input of type %s not supported' % (
                    type(pareto), type(v)))
            vec = oalib.mvalue_t_long(list(v))
            vector_pareto.push_back(vec)
    elif isinstance(pareto, oalib.ParetoMultiDoubleLong):
        vector_pareto = oalib.vector_mvalue_t_double()
        for v in values:
            if isinstance(v, (int, float)):
                # convert to list type
                v = [float(v)]
            if not isinstance(v, (list, tuple)):
                raise Exception('creating Pareto element for Pareto object of type %s and input of type %s not supported' % (
                    type(pareto), type(v)))

            vec = oalib.mvalue_t_double(list(v))
            vector_pareto.push_back(vec)
    elif isinstance(pareto, oalib.ParetoDoubleLong):
        if not isinstance(values, (list, tuple, np.ndarray)):
            raise Exception('cannot handle input of type %s' % (tuple(values), ))
        vector_pareto = values
    else:
        raise Exception('creating Pareto element for Pareto object of type %s and input of type %s not supported' % (
            type(pareto), type(values)))
    return vector_pareto
