pde-calculator
==============

A Finite Difference Calculator for solving Hyperbolic PDEs.

This code numerically solves hyperbolic PDEs of the form:  

    Dt[u(x, t)] + a(t, x) * Dx[u(x, t)] = F(x, t) where 
    
    Dt[] and Dx[] are the differential operators for t and x

The solutions are animated in a window. Data is saved to text files when 's' is pressed.

Schemes
-------

The following finite difference schemes are implemented:

* Forward-Time Forward-Space
* Forward-Time Backward-Space
* Forward-Time Central-Space
* Lax Fredrichs
* Leapfrog
* Equilibrium
* Lax Wendroff

