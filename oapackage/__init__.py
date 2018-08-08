""" Orthogonal Array package

The Orthogonal Array package is a pacakge to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""

import oalib
import oapackage.Doptim
import oapackage.tests

oalib.setloglevel(oalib.SYSTEM)
oalib.log_print(-oalib.SYSTEM, '')

try:
    from oalib import *
except:
    # fix for RTD
    pass
from . oahelper import *
from . Doptim import *
from . import scanf

import numpy as np

__all__ = ['oahelper']
__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
__doc__ = __description__ + " <" + __uri__ + ">"

__version__ = oalib.version()
