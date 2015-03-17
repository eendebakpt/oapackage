# OApackage
# pieter.eendebak@gmail.com

__all__ = ['oahelper']

import oalib
oalib.setloglevel(oalib.SYSTEM)
oalib.log_print(-oalib.SYSTEM, '')
from oalib import *
from . oahelper import *
from . Doptim import *

#%%

def unittest():
  print('oapackage: oalib version %s' % oalib.version() )
  al=oalib.array_link()
  ii=0
  al=oalib.exampleArray(ii, 0)
  Deff=al.Defficiency()
  print('## oapackage test: example array %d: Deff %.3f' % (ii, Deff))

  