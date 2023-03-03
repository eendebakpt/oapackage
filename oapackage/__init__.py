import oapackage.conference
import oapackage.Doptim
import oapackage.oahelper
from oapackage.Doptim import (  # noqa
    Doptimize,
    array2Dtable,
    calcScore,
    filterPareto,
    generateDscatter,
    optimDeffPython,
    scoreDn,
    selectDn,
)
from oapackage.oahelper import (  # noqa
    DefficiencyBound,
    array2html,
    array2latex,
    arrayfile_generator,
    checkArrayFile,
    checkFiles,
    checkFilesOA,
    checkOAfile,
    choose,
    compressOAfile,
    create_pareto_element,
    createPareto,
    deprecated,
    designStandardError,
    enlargelims,
    extendSingleArray,
    fac,
    fileinput,
    finddirectories,
    findfiles,
    findfilesR,
    floatformat,
    formatC,
    getArrayFile,
    guitest_monitorSizes,
    gwlp2str,
    helmert_contrasts,
    isstr,
    joinArrayLists,
    jseq,
    makearraylink,
    mkdirc,
    monitorSizes,
    nArrayFile,
    niceplot,
    oainfo,
    oaIsBinary,
    oalib,
    parseProcessingTime,
    plot2Dline,
    randomizearrayfile,
    runcommand,
    runExtend,
    safemax,
    safemin,
    selectArrays,
    selectArraysInFile,
    selectJ,
    selectParetoArrays,
    series2htmlstr,
    setWindowRectangle,
    sortcols,
    sortrows,
    static_var,
    strftime,
    test_selectJ,
    tilefigs,
    timeString,
    tprint,
    webbrowser,
    write_text_arrayfile,
)

""" Orthogonal Array package

The Orthogonal Array package is a package to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""
from oalib import *  # noqa # legacy structure #type: ignore
from oalib import (  # noqa #type: ignore
    ParetoDoubleLong,
    array_link,
    arraydata_t,
    arraylink2arraydata,
    exampleArray,
    reduceGraphNauty,
    arrayfile_t,
    reduceOAnauty,
    transformGraphMatrix,
)

oapackage.oalib.setloglevel(oapackage.oalib.SYSTEM)
oapackage.oalib.log_print(-oapackage.oalib.SYSTEM, "")
__version__ = oapackage.oalib.version()

__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
