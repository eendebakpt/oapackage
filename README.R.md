Orthogonal Array Package
========================

An example file for using the Orthogonal Array packages is contained in scripts/example.R

Installation
------------

To build the package from the source code use:

> R CMD build oapackage

Installation from a package file can be done with

> R CMD INSTALL [PACKAGEFILE]


Development
-----------

> install.packages('devtools')
> library(devtools)
> load_all('oapackage')
> document("oapackage") 