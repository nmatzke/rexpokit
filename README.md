# rexpokit

Wraps some of the matrix exponentiation utilities from
[EXPOKIT](http://www.maths.uq.edu.au/expokit/), 
a Fortran 77 library that is widely recommended for matrix exponentiation 
(Sidje RB, 1998. "Expokit: A Software Package for Computing Matrix 
Exponentials." ACM Trans. Math. Softw. 24(1): 130-156). EXPOKIT includes 
functions for exponentiating both small, dense matrices, and large, sparse 
matrices (in sparse matrices, most of the cells have value 0). 

Rapid matrix exponentiation is useful in phylogenetics when we have a large 
number of states (as we do when we are inferring the history of transitions 
between the possible geographic ranges of a species), but is probably 
useful in other ways as well.

Build status on Travis-CI:

[![Build Status](https://travis-ci.org/nmatzke/rexpokit.svg?branch=master)](https://travis-ci.org/nmatzke/rexpokit)
