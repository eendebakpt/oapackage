#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Example script for Python interface to Orthogonal Array code.

mathescape, numbersep=5pt, label=Python session, gobble=0, frame=lines, framesep=1.5mm, bgcolor=bg

> pygmentize -P "bgcolor=.95,.95,.95" -P "style=colorful" -f html example_python_testing.py
> pygmentize -S colorful -f html > style.css; cp style.css ~/misc/homepage/oapage/


@author: Pieter Eendebak
"""


#import oalib
import oapackage

print('oalib version: %s' % oapackage.version() )

al=oapackage.exampleArray(0)
al.showarray()
print('D-efficiency %f, rank %d' % (al.Defficiency(), al.rank()) )
print('Generalized wordlength pattern: %s' % str(al.GWLP()))


al=oapackage.exampleArray(1)
print('Generalized wordlength pattern: %s' % str(al.GWLP()))


