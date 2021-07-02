from . import scanf
from .Doptim import *  # type: ignore
from .oahelper import *

""" Orthogonal Array package

The Orthogonal Array package is a package to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""

import oalib
import oapackage.Doptim
from oalib import *
from oalib import (ParetoDoubleLong, array_link, arraydata_t, arrayfile_t,
                   arraylink2arraydata, exampleArray, reduceGraphNauty,
                   reduceOAnauty, transformGraphMatrix)

from . import conference

oapackage.oalib.setloglevel(oapackage.oalib.SYSTEM)
oapackage.oalib.log_print(-oapackage.oalib.SYSTEM, '')
__version__ = oapackage.oalib.version()

__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
