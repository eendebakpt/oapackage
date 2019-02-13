---
title: 'OApackage: A Python package for generation and analysis of orthogonal arrays, optimal designs and conference designs'
tags:
  - orthogonal arrays
  - optimal designs
  - conference matrices
  - conference designs
  - design of experiments
authors:
  - name: Pieter Thijs Eendebak
    orcid: 0000-0001-7018-1124
    affiliation: "1, 2"
  - name: Alan Roberto Vazquez
    orcid: 0000-0002-3658-0911
    affiliation: "1, 3"
affiliations:
 - name: Department of Engineering Management, University of Antwerp, Belgium
   index: 1
 - name: Netherlands Organisation for Applied Scientific Research (TNO), P.O. Box 155, 2600 AD Delft, The Netherlands
   index: 2
 - name: Department of Biosystems, KU Leuven, Leuven, Belgium
   index: 3
date: 11 November 2018
bibliography: paper.bib
---

# Summary

Orthogonal arrays, optimal designs and conference designs are important tools for the design of
experiments [@Elster1995], [@hedayat2012orthogonal], [@Wu2009]. The OApackage (Orthogonal Array package) contains functionality 
to generate and analyse these types of designs. More specifically, the OApackage allows 
the user to:

* Efficiently generate orthogonal arrays, optimal designs and conference designs
* Reduce the designs to their normal form and perform isomorphism testing 
* Calculate a wide variety of statistical properties of the designs

The data analysis of the experiments conducted using the generated designs is left to 
existing statistical software such as R [@Rpackage] and JMP [@wiki:JMP].

To generate orthogonal arrays and conference designs, the OApackage uses an exhaustive 
generation procedure with isomorphism pruning [@Eendebak2009], [@Schoen2019]. To generate 
optimal designs, the package uses a flexible optimality criterion and a coordinate-exchange 
optimization algorithm [@Eendebak2015].

The reduction of the designs to their normal form is done by either reduction to a minimal form 
(such as lexicographically minimal in columns or delete-one-factor projection normal form [@EendebakDOF])
or reduction using graph algorithms. For designs with a specified isomorphism group,
the OApackage provides a generic interface to the graph reduction algorithms that effectively perform isomorphism 
testing and reduction to normal form.

The OApackage evaluates the orthogonal arrays, optimal designs and conference designs using well-known statistical criteria. For instance,
the package can calculate the generalized wordlength pattern and confounding frequency vector [@Tang1999], which are based
on the J-characteristics [@Deng1999],
and the number of degrees of freedom available for estimating selected factors' effects.
The package can also calculate the $F_4$ vector of
conference designs [@Schoen2019] and the D-efficiency of optimal designs [@Goos2011]. 

The OApackage consists of a C++ library with a Python interface generated
by SWIG. The source code is available at https://github.com/eendebakpt/oapackage. Examples for both 
generation and analysis of designs are available in the OApackage documentation [@OAdocumentation].
The Orthogonal Array package website [@EendebakOAwebsite] contains a large collection of orthogonal arrays, 
optimal designs and conference designs generated with the package.

# Acknowledgements

We acknowledge useful discussions with Eric Schoen during the development of this project.

# References
