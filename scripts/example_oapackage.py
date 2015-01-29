#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Example script for Python interface to Orthogonal Array code.

@author: Pieter Eendebak
"""

import oapackage
print('oalib version: %s' % oapackage.version() )

al=oapackage.exampleArray(0)
al.showarray()
print('D-efficiency %f, rank %d' % (al.Defficiency(), al.rank()) )
print('Generalized wordlength pattern: %s' % str(al.GWLP()))

al=oapackage.exampleArray(1)
print('Generalized wordlength pattern: %s' % str(al.GWLP()))


