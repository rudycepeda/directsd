function P = findzero ( P, tol )
%FINDZERO  Zeroing very small elements of a polynomial.
%       
%     Q = FINDZERO ( P, TOL )
%
%   Inputs:
%     P - polynomial
%     TOL - tolerance (default eps)
%
%   Outputs:
%     Q - reduced polynomial
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        if ~exist('tol', 'var'), tol = eps; end;
        tol = tol * norm(P);
        if norm(P) < 1.e-10, tol = 1.e-10; end;
        if tol < eps, tol = eps; end;

        ind = find(abs(P) < tol);
        if ~isempty(ind)
           P(ind) = 0;
           P = striplz(P);
        end;
        
%------- End of FINDZERO.M --------- KYuP ----------
