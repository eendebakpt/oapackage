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

def autodoctest():
    """ Test the module using autodoc
    Example:   
      >>> import oapackage
      >>> arrayclass=oapackage.arraydata_t(2, 40, 0, 7)
      >>> print(arrayclass)
      arrayclass: N 40, k 7, strength 0, s {2,2,2,2,2,2,2}, order 0

    """
    return

def unittest():
  print('oapackage: oalib version %s' % oalib.version() )
  al=oalib.array_link()
  ii=0
  al=oalib.exampleArray(ii, 0)
  Deff=al.Defficiency()
  print('## oapackage test: example array %d: Deff %.3f' % (ii, Deff))

if __name__ == "__main__":
    """ Dummy main for oapackage """
    import doctest
    doctest.testmod()    
  