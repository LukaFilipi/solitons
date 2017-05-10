# Solitons
Simulation of soliton dynamics

C++ simulation of propagation of solitons according to the Korteweg-de Vries (KdeV) equation.

A centre-difference approximation is used for space and a 4th order Runge Kutta method is used for time discretisation.

The theory is described in detail in the file solitons.pdf

A range of different initial solutions were propagated through the KdeV equation to study the behaviour of soliton dynamics, including:
1. Normal Mode - a single soliton
2. Collisions - 2 normal mode solutions starting at different places
3. Wave breaking - a solution that is not a normal mode
4. Shock waves - the dispersion term is removed from the KdeV equation

The effects of changing some of the equation parameters on the stability of the numerical solution are also examined.

The directories inside the "Results" directory contain the plotted graphs for each of these scenarios with different parameters. The different parameters used for each of the plots are indicated in the filenames.
