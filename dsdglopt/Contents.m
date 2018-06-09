% DirectSD Toolbox - Global optimization
% Version 3.0  01-May-2006
%
%   admproj     - Proection of poles onto an admissible region.
%   arandsearch - Accelerated random search optimization.
%   banana      - Rosenbrock's "banana function" f(x,y)= 100*(y-x^2)^2+(1-x)^2.
%   bin2val     - Transform a binary string value to a real number.
%   coord2hilb  - Map an ND-vector to 1D-value in [0,1] using Hilbert curve.
%   cp2par      - Map a polynomial to an array of parameters in [0,1]. 
%   crandsearch - Controlled random search (CRS) optimization.
%   f_sdh2p     - Target function for nonlinear H2-optimization.
%   f_sdl2p     - Target function for nonlinear L2-optimization.
%   go_par2k    - Map array of parameters in [0,1] to the controller.
%   go_sdh2p    - Target function for modal-based H2-optimization.
%   go_sdl2p    - Target function for modal-based L2-optimization.
%   guesspoles  - Select N poles as an initial guess.
%   hilb2coord  - Map an 1D-value in [0,1] to ND-vector using Hilbert curve.
%   infglob     - Global optimization using information algorithm.
%   infglobc    - Constrained global optimization using information algorithm. 
%   k2ksi       - Find a polynomial parameter \xi for a controller.
%   neldermead  - Simplex minimization method by Nelder and Mead.
%   optglob     - Global optimization using information algorithm.
%   optglobc    - Constrained global optimization using information algorithm.
%   par2cp      - Map array of parameters in [0,1] to the characteristic polynomial.
%   r1range     - Find the range of r1.
%   r2range     - Find the range of r2.
%   randbeta    - Beta-distributed random numbers.
%   randgamma   - Gamma-distributed random numbers.
%   randsearch  - Random search optimization.
%   sa_testfun  - Complex function with many local minima.
%   sasimplex   - Simulated annealing using Nelder-Mead method. 
%   sector      - Computes stability degree abs relative oscillation.
%   simanneal   - Direct search using simulated annealing technique.
%   u2range     - Linear mapping from [0,1] onto [min,max].
%   uniproj     - Project values onto [0,1].
%   updateopt   - Update options structure.
%   val2bin     - Map a real value to a binary one.

%   Copyright (c) 1995-2006 by Konstantin Polyakov
%   $Revision: 3.00 $  $Date: 07-Apr-2006 $
