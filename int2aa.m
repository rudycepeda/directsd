function I = int2aa ( A, T )
%INT2AA Compute double integral involving matrix exponential.
%
%     I = INT2AA ( A, T )
%
%   Inputs:
%     A - square matrix
%     T - upper integration limit
%   Outputs:
%     I - value of the integral:
%           I = int[0..T] exp(A*t) int[0..t] exp(-A*v) dv dt
%
%   See also INTAXA, INTEAT.

% References:
% [1] C.F. Van Loan, Computing integrals involving the matrix
%     exponentials, IEEE Trans. Autom. Contr., vol. AC-23, no. 3,
%     June 1978, pp. 395-404.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 03-May-2006 $
%------------------------------------------------------

        n = size(A, 1);
        In = eye(n);
        On = zeros(n);
        C = [-A In On
              On On In
              On On -A];
        E = expm(C*T);
        I = E(1:n,1:n) \ E(1:n,2*n+1:end);

%------- End of INT2AA.M --------- KYuP ----------
