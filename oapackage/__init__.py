""" Orthogonal Array package

The Orthogonal Array package is a pacakge to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""

import oapackage.oalib
import oapackage.Doptim


from . oalib import arraydata_t, array_link, exampleArray, ParetoDoubleLong, reduceOAnauty, arraylink2arraydata, reduceGraphNauty, transformGraphMatrix

try:
    from oalib import *
except:
    # fix for RTD
    pass

from . oahelper import *
from . Doptim import *
from . import scanf

import numpy as np

oalib.setloglevel(oalib.SYSTEM)
oalib.log_print(-oalib.SYSTEM, '')

#__all__ = ['oahelper']
__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
__doc__ = __description__ + " <" + __uri__ + ">"

__version__ = oalib.version()
