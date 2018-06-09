function [A,nz] = delzero ( A, tol )
%DELZERO Delete trailing zero members of A.
%
%     [B,nz] = DELZERO ( A, TOL )
%
%   Inputs:
%     A - a polynomial
%     TOL - tolerance (default is eps)
%
%   Outputs:
%     B  - reduced polynomial
%     nz - the number of zeros deleted
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('tol','var'), tol = eps; end;
        nz = 0;
        if isempty(A), return; end;
%------------------------------------------------------
%       Polynomial case
%------------------------------------------------------
        if isa(A,'poln')
           z = A.z; 
           ind = find(abs(z) < tol);
           nz = length(ind);
           if nz == 0, return; end;
           z(ind) = [];
           coef = A.coef(1:end-nz);
           A = change(A, A.k, z, coef);
           return;
        end;
%------------------------------------------------------
%       Numerical case
%------------------------------------------------------
        A = striplz(A);
        i = size ( A, 2 );
        while i > 1
           if norm(A(i)) > tol, break; end;
           A = A(1:i-1);
           nz = nz + 1;
           i = i - 1;
        end;

%------- End of DELZERO.M --------- KYuP ----------
