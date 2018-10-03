# rexpokit

Wraps some of the matrix exponentiation utilities from
[EXPOKIT](http://www.maths.uq.edu.au/expokit/), 
a Fortran 77 library that is widely recommended for matrix exponentiation 
(Sidje RB, 1998. "Expokit: A Software Package for Computing Matrix 
Exponentials." *ACM Trans. Math. Softw.* 24(1): 130-156). EXPOKIT includes 
functions for exponentiating both small, dense matrices, and large, sparse 
matrices (in sparse matrices, most of the cells have value 0). 

Rapid matrix exponentiation is useful in phylogenetics when we have a large 
number of states (as we do when we are inferring the history of transitions 
between the possible geographic ranges of a species), but is probably 
useful in other ways as well.

**Build status** on Travis-CI: [![Build Status](https://travis-ci.org/nmatzke/rexpokit.svg?branch=master)](https://travis-ci.org/nmatzke/rexpokit)

**NOTE:** As of 2018-10-03, a new version of rexpokit, 0.26.6, has been accepted on CRAN. This version fixes warnings due to upgrades in CRAN's FORTRAN compilers, and more importantly, removes some dependencies previously needed by cladoRcpp/BioGeoBEARS. CRAN version is here, binaries should be available in a few days: [https://cran.r-project.org/package=rexpokit](https://cran.r-project.org/package=rexpokit)

**Release v0.26.6** registered on Zenodo: [![DOI](https://zenodo.org/badge/17001945.svg)](https://zenodo.org/badge/latestdoi/17001945)

**Zenodo link for release:** https://zenodo.org/badge/latestdoi/17001945

**Zenodo DOI for release:** http://dx.doi.org/10.5281/zenodo.1442889
