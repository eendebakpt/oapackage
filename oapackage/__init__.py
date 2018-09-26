""" Orthogonal Array package

The Orthogonal Array package is a package to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""

import oalib
import oapackage.Doptim
from oalib import *
from oalib import arraydata_t, array_link, exampleArray, ParetoDoubleLong, reduceOAnauty, arraylink2arraydata, reduceGraphNauty, transformGraphMatrix

oapackage.oalib.setloglevel(oapackage.oalib.SYSTEM)
oapackage.oalib.log_print(-oapackage.oalib.SYSTEM, '')
__version__ = oapackage.oalib.version()

from . oahelper import *
from . Doptim import *
from . import scanf


__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
#__doc__ = __description__ + " <" + __uri__ + ">"
