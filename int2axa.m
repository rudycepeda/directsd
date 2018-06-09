function I = int2axa ( A, X, T )
%INT2AXA Compute double integral involving matrix exponential.
%
%     I = INT2AXA ( A, X, T )
%
%   Inputs:
%     A, X - square matrices
%     T - upper integration limit
%   Outputs:
%     I - value of the integral:
%           I = int[0..T] int[0..t] exp(A*v)*X*exp(A'*v) dv dt
%
%   See also INT2AA, INTAXA, INTEAT.

% References:
% [1] C.F. Van Loan, Computing integrals involving the matrix
%     exponentials, IEEE Trans. Autom. Contr., vol. AC-23, no. 3,
%     June 1978, pp. 395-404.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 26-Sep-2006 $
%------------------------------------------------------

        n = size(A, 1);
        In = eye(n);
        On = zeros(n);
        C = [-A In On
              On -A X
              On On A'];
        E = expm(C*T);
        I = E(1:n,1:n) \ E(1:n,2*n+1:end);

%------- End of INT2AXA.M --------- KYuP ----------
