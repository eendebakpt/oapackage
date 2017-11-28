Orthogonal Array Package
========================

An example file for using the Orthogonal Array packages is contained in tests/example_Doptimize.R

Installation
------------

To install the package from CRAN use (from within R):

>  install.packages('oapackage')

To build the package from the source code use:

> R CMD build oapackage

Installation from a package file can be done with

> R CMD INSTALL [PACKAGEFILE]

Make sure yoy have rcppeigen installed.

Development
-----------

> install.packages('devtools')
> library(devtools)
> load_all('oapackage')
> document("oapackage") 
