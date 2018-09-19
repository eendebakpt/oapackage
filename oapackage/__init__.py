""" Orthogonal Array package

The Orthogonal Array package is a pacakge to generate and analyse orthogonal
arrays, optimal designs and conference matrices. For more information see

http://github.com/eendebakpt/oapackage

"""

import os
_rtd = os.environ.get('READTHEDOCS', False)

if _rtd:
    # fix for RTD
    oalib = None
    __version__ = 'RTD'
else:
    import oapackage.oalib as oalib
    import oapackage.Doptim
    from . oalib import *

    oapackage.oalib.setloglevel(oapackage.oalib.SYSTEM)
    oapackage.oalib.log_print(-oapackage.oalib.SYSTEM, '')
    __version__ = oapackage.oalib.version()

    from . oalib import arraydata_t, array_link, exampleArray, ParetoDoubleLong, reduceOAnauty, arraylink2arraydata, reduceGraphNauty, transformGraphMatrix
    
    from . oahelper import *
    from . Doptim import *
    from . import scanf


__description__ = "Orthogonal Array package"
__uri__ = "http://www.pietereendebak.nl/oapackage/index.html"
#__doc__ = __description__ + " <" + __uri__ + ">"

