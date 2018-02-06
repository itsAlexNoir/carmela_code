# carmela_code
## A python library for solving the one-electron one-dimensional TDSE.

This library solves the time-dependent Schoedinger equation (TDSE) in
a one-electron one-dimensional toy model. The purpose of the library
was to easily play with a TDSE for fast development. The library and
modules therein are written in python 3.

---

The library can be imported with the command::

`from carmela_code import *`

The modules included are the following:

   :mod:`absorber`
   :mod:`axes`
   :mod:`constants`
   :mod:`fdrule`
   :mod:`flux`
   :mod:`graph`
   :mod:`hamiltonian`
   :mod:`input_reader`
   :mod:`laser`
   :mod:`momentum`
   :mod:`wavefunction`

In order to run the library, the following packages are required: 
numpy, matplotlib, palettable, fileinput and sys. 
Except palettable the rest of the packages are available in most of the current python distributions.
