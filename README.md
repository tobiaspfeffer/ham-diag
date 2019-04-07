This repository provides an algorithm for the solution of the non-linear
integral equation of the form

f(x) = c(x) + int_0^1 d^Dy K(x,y) f(y)^2

for given functions c, K, and x,y \in \mathcal{R}^D is a D dimensional vector.
The detailed algorithm is described in Chapter 5.4 of my PhD thesis available
at https://edoc.ub.uni-muenchen.de/23644/.

The novelty of the algorithm is that it introduces the first and only method to
avoid the curse of dimensionality in the solution of these kind of integral
equations, i.e., the scaling of the computational time is independent of the
dimensionality D. Whereas, e.g., for grid based methods, the computational
time scales generically exponential with the dimensionality,
i.e., \mathcal{O}(\exp(D)), making it impossible for state-of-the-art methods
to go beyond D>3. This novel algorithm provides a stochastic solution through
a Markov Chain Monte Carlo process summing up a homotopy series for the solution
of the integral equation. Therefore, as generically for Monte Carlo methods, the
algorithm scales with the number of stochastically independent samples N
collected in the Markov chain process, \mathcal{O}(N^1/2).
For further details consult Chapter 5.4 of my PhD thesis available at 
https://edoc.ub.uni-muenchen.de/23644/.

General Usage
-------------

License
-------

Copyright © 2019  Tobias Pfeffer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the
file [LICENSE.txt](LICENSE.txt).