import oapackage.Doptim

from .Doptim import *  # noqa
from .oahelper import *  # noqa

""" Orthogonal Array package

The Orthogonal Array package is a package to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""

from oalib import *  # noqa # legacy structure
from oalib import arrayfile_t  # noqa
from oalib import reduceGraphNauty  # noqa
from oalib import (  # noqa
    ParetoDoubleLong,
    array_link,
    arraydata_t,
    arraylink2arraydata,
    exampleArray,
    reduceOAnauty,
    transformGraphMatrix,
)

oapackage.oalib.setloglevel(oapackage.oalib.SYSTEM)
oapackage.oalib.log_print(-oapackage.oalib.SYSTEM, "")
__version__ = oapackage.oalib.version()

__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
