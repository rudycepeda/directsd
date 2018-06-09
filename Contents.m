% DirectSD Toolbox
% Version 3.0  01-May-2006
%
% Polynomials.
%   s           - Construct elementary polynomial: s + 0.
%   q           - Construct elementary polynomial: q + 0.
%   p           - Construct elementary polynomial: p + 0.
%   d           - Construct elementary polynomial: d + 0.
%   z           - Construct elementary polynomial: z + 0.
%   compat      - Make compatible polynomials with the same variable.
%   deg         - Degree of a polynomial.
%   roots       - Roots of a polynomial.
%   polyval     - Value of a polynomial.
%   derive      - Derivative of a polynomial.
%   polyder     - Value of a polynomial derivative.
%   isct        - Returns 1 for a continuous-time polynomial and 0 otherwise.
%   isdt        - Returns 1 for a discrete-time polynomial and 0 otherwise.
%   gcd         - Monic greatest common divisor of two or more polynomials.
%   tf2nd       - Extract numerator and denominator polynomials.
%   coprime     - Extract greatest common divisor of two polynomials.
%   triple      - Extract greatest common divisor of there polynomials.
%   factor      - Stable-unstable factorization of a polynomial.
%   striplz     - Strip leading zeros of a polynomial.
%   sfactor     - Spectral factorization for polynomials and transfer functions.
%   sfactfft    - Polynomial spectral factorization using FFT.
%   recip       - Reciprocal polynomial.
%   c2z         - Discretization into the z-plane.
%   c2d         - Discretization into the d-plane (d = 1/z).
%   delzero     - Extarct zeros at the origin.
%   zpk         - Transform to a zpk-model.
%   vec         - Vectorize coefficients of a (quasi)polynomial.
%
% Polynomial equations.
%   dioph       - Simple Diophantine polynomial equation.
%   dioph2      - Solves a special polynomial equations X*A + X~*B + Y*C = 0.
%   diophsys    - Solves a system of two Diophantine polynomial equations.
%   diophsys2   - Solves a system of two special polynomial equations.
%
% Discrete trasformations.
%   ztrm        - Modified Z-transform for a transfer matrix.
%   dtfm        - Discrete Laplace transform with a hold.
%   dtfm2       - Special symmetric discrete Laplace transform with a hold.
%
% Analysis of sampled-data systems.
%   charpol     - Characteristic polynomial of sampled-data system.
%   sdmargin    - Stability margin for sampled-data system.
%   quaderr     - Squared integral error for sampled-data system.
%   sdl2err     - Integral quadratic error for sampled-data systems.
%   sdh2norm    - H2-norm of sampled-data systems.
%   sd2doferr   - Integral quadratic error in 2-DOF sampled-data system.
%   dinfnorm    - Hinf-norm for discrete-time systems.
%   dahinorm    - Associated Hinf-norm for discrete-time system.
%   sdhinorm    - Hinf-norm for sampled-data systems.
%   sdahinorm   - Associated Hinf-norm for sampled-data system.
%   sdtrhinferr - Associated Hinf-cost for sampled-data tracking systems.
%   sdnorm      - Various norms of sampled-data systems.
%   sdsim       - Simulation of sampled-data systems.
%
% Polynomial design.
%   ch2         - H2-optimal continuous-time controller.
%   dhinf       - Polynomial Hinf-optimal discrete-time system design.
%   polquad     - Polynomial minimizition of a quadratic functional.
%   whquad      - Wiener-Hopf minimization of frequency-domain quadratic functionals.
%   sdh2        - H2-optimal controller for sampled-data systems.
%   sdh2hinf    - Mixed H2/AHinf-optimization of sampled-data systems.
%   sdl2        - L2-optimal controller for sampled-data systems.
%   sd2dof      - Optimal feedforward controller for 2-DOF sampled-data systems.
%   split2dof   - Stable realization of 2-DOF sampled-data systems.
%   sdahinf     - AHinf-optimal controller for sampled-data systems.
%   polhinf     - Polynomial solution to an Hinf-problem.
%   sdtrhinf    - AHinf-optimal controller for sampled-data tracking systems.
%   polopth2    - Solution to polynomial H2-minimization problem. 
%   modsdh2     - H2-optimal reduced-order modal digital controller.
%   modsdl2     - L2-optimal reduced-order modal digital controller.
%   psigain     - DC-gain of Psi-parameter for tracking system.
%
% State-apace design.
%   ssquad      - State-space minimization of frequency-domain quadratic functionals.
%   h2reg       - H2-optimal controller for LTI systems.
%   hinfreg     - Hinf-optimal controller for LTI systems.
%   sdfast      - Fast discretization of sampled-data systems.
%   sdgh2mod    - Equivalent discrete model in generalized H2-problem. 
%   sdh2simple  - Equivalent discrete model in simple H2-problem. 
%   sdh2reg     - H2-optimal controller for sampled-data systems.
%   sdhinfreg   - Hinf-optimal controller for sampled-data systems. 
%
% Miscellaneous functions.
%   bilintr     - Bilinear transformation for SISO system.
%   improper    - Separate improper (polynomial) part of a rational matrix.
%   sector      - Find degree of stability and oscillation for closed-loop poles.
%   sumzpk      - Reliable summation of zpk-models with common poles.
%   sfactor     - Spectral factorization for polynomials and transfer functions.
%   separss     - Proper separation (state-space technique).
%   separtf     - Proper separation (polynomial equations technique).
%
% Linear algebra.
%   Type "help dsdlinalg" for a list of linear algebra functions.
%
% Global optimization.
%   Type "help dsdglopt" for a list of global optimization functions.
%
% Demonstrations.
%   Type "demo" or "help dsddemos" for a list of available demos.

%   Copyright (c) 1995-2006 by Konstantin Polyakov
%   $Revision: 3.00 $  $Date: 07-Apr-2006 $
