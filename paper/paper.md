---
title: 'OApackage: A Python package for generation and analysis of orthogonal arrays and conference designs'
tags:
  - Python
  - ortogonal arrays
  - conference matrices
  - design of experiments
authors:
  - name: Pieter Thijs Eendebak
    orcid: 0000-0001-7018-1124
    affiliation: "1, 2"
  - name: Alan Roberto Vazquez
    orcid: 0000-0002-3658-0911
    affiliation: "1, 3"
affiliations:
 - name: University of Antwerp
   index: 1
 - name: TNO
   index: 2
 - name: University of Leuven
   index: 3
date: 27 September 2018
bibliography: paper.bib
---

# Summary

Orthogonal arrays, optimal designs and conference designs are important tools for the design of
experiments [@Elster1995], [@Hedayat1999], [@Wu2009]. The Orthogonal Array package contains functionality 
to generate and analyse these types of designs. More specifically, the Ortogonal Array package allows 
the user to:

* Efficiently generate orthogonal arrays, optimal designs and conference designs
* Reduce the designs to their normal form and perform isomophism testing 
* Calculate a wide variety of statistical properties of the designs

To generate orthogonal arrays and conference designs, the Orthogonal Array package uses an exhaustive 
generation procedure with ismorphism pruning [@Eendebak2009], [@IsomorphismPaper]. To generate 
optimal designs, the package uses a flexible optimality criterion and a coordinate-exchange 
optimization algorithm [@Eendebak2015].

The reduction of the designs to their normal form is done by either reduction to a minimal form 
(such as lexicographic minimal in columns or delete-one-factor projection normal form [@EendebakDOF])
or reduction using graph algorithms. For designs with a specified isomorphic group, the Orthogonal 
Array package provides a generic interface to the graph reduction algorithms to perform ismophism 
testing and reduction to normal form effectively.

The Orthogonal Array package evaluates the designs using sensible statistical criteria. For instance,
the package calculates the generalized wordlength pattern [@Tang1999], the J_k-characteristics [@Deng1999] 
and the number of degrees of freedom available for estimating selected factors' effects.

The Orthogonal Array package consists of a C++ library with an user-friendly Python interface generated
by SWIG. The source code is available at https://github.com/eendebakpt/oapackage. Examples for both 
generation and analysis of designs are available in the OApackage documentation [@OAdocumentation].

The Orthogonal Array package website [@EendebakOA] contains a large collection of orthogonal arrays, 
optimal designs and conference designs. An alternative collection of orthogonal arrays can be found in 
the website of Neil Sloane [@Sloanewebsite]. Finally, the analysis of data from the designs is left to 
exising packages such as R [@Rpackage] and JMP [@wiki:JMP].

# Acknowledgements

We acknowledge useful discussions with Eric Schoen during the development of this project.

# References
