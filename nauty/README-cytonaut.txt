Cytonaut - A Python Interface for Nauty
=======================================
This package is a quick-and-dirty Python interface for Brendan McKay
and Adolfo Piperno's "nauty" tool. Nauty must be downloaded separately
from https://cs.anu.edu.au/people/Brendan.McKay/nauty/ and is subject
to its own terms of use.

Cytonaut also requires the NetworkX library, available at
https://networkx.github.io/ , and expects input to be in the form of
NetworkX graphs.

Cytonaut is written in Cython; while compiling it out-of-the-box
does not require Cython (it ships with Cython-generated code), you
will need Cython if you wish to change "cytonaut.pyx" and recompile.
You can find Cython at http://cython.org/ .

How to build
============
In order to build cytonaut, the cytonaut directory must have a copy of
the nauty source code in a subdirectory named "nauty". The easiest way
to do this is to create a symlink in the cytonaut directory pointing to
the location of the nauty source code.

Before building cytonaut, you must build nauty, so that its object files
are still available.

To build Cytonaut from the prepackaged generated code, run "make" in its
directory. To build cytonaut from its Cython source (requires Cython),
instead run "make cython".

After building Cytonaut, you may install it as a site package by running
"make install".

Using Cytonaut
==============
The file "example.py" contains a simple example of Cytonaut usage.

In general, given a NetworkX graph G, you will want to construct a
CytoGraph object, cytonaut.CytoGraph(G). The CytoGraph object supports
calls to nauty. Note that a CytoGraph object created this way will
*NOT* track any changes to the original graph, and should effectively
be considered an immutable copy.

If you wish to specify a node coloring of G which the automorphisms
and canonical form must respect, you should use the node attribute
'color'.  An example of this is shown in example.py . Nodes with no
color specified are all given the color 0 by default.

To run nauty on a CytoGraph object CG, call "CG.do_nauty()". You
may optionally specify the keyword argument "aut", which, if True,
instructs nauty to compute and store the automorphism group of CG.

After calling do_nauty, you may access the various computed properties
of the graph through the methods "canon_string", "canon_translation",
"orbits", and "automorphisms", as demonstrated in example.py .

Known Issues
============
Cytonaut may fail to work on graphs with more than 32 or 64 vertices
(depending on your system's word size).

Licensing
=========
This project has been released into the public domain via a CC0 license.
See COPYING.txt for legalese.

Disclaimer
==========
I make no guarantees of correctness (or prettiness) of this software
and I take no responsibility for any damage it may cause.
