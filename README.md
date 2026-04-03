# 1D_Linear-Hyp
Solver for Simple Linear Hyperbolic Equation with constant velocity. The program contains the following algorithms - 
* Upwind Scheme
* Lax-Wendroff Scheme
* MacCormack Scheme

The results have been worked out for the initial condition ```exp(-100(x-0.5)^2)``` with periodic boundary conditions and the results for the three methods compared after a single cycle by evaluating the absolute errors against the initial values. The results have been plotted in Python using ```matplotlib```.

The code has three sections, with the variables (Datastruct.cpp) and functions (Functions.cpp) being define in separate files, while the main function (Combined_Schemes.cpp) containing only the functions. The parameters are extracted from ```Read_File.txt```, which can be changed as per requirement.

Any queries can be forwarded to **shounakd@jncasr.ac.in**

### References - 
* Chattot, Jean-Jacques. Computational Aerodynamics and Fluid Dynamics: An Introduction. Scientific Computation. Springer, 2002. https://doi.org/10.1007/978-3-662-05064-4.
