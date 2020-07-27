# Notice:
This is a modification of Kaa for experimenting with other optimization techniques besides Bernstein polynomials. Mostly trying to find ways to extend beyond polynomial dynamics.
This version utilizes a python wrapper over [Kodiak](https://github.com/nasa/Kodiak). Many thanks to [Stanley Bak ](https://github.com/stanleybak) for providing the original wrapper over Kodiak.

# Kaa
Kaa is a tool for reachability analysis of polynomial dynamical systems using parallelotope bundles.
It is a rewrite of the tool Sapo introduced by Dreossi, Dang, Piazza [(paper)](https://dl.acm.org/doi/abs/10.1145/2883817.2883838)

# Dependencies
Kaa relies only on the following python3 packages:

- [sympy](https://pypi.org/project/sympy/)
- [scipy](https://pypi.org/project/scipy/)
- [numpy](https://pypi.org/project/numpy/)
- [swiglpk](https://pypi.org/project/swiglpk/)

All of which can be installed through pip or through the package's corresponding page on PyPI.

# Running Sample Models and Examples.

The Juypter notebook named kaa-intro provides basic examples of the usage of kaa. The notebook contains all of the relevant code necessary to begin plotting the reachable sets and phase plots.

# Contents:
- [Summary of Files](md/explan.md)
