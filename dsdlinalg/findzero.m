function A = findzero ( A, tol )
%FINDZERO Make zero all negligible elements of a matrix.
%
%     A0 = FINDZERO ( A )
%     A0 = FINDZERO ( A, TOL )
%
%   Inputs: 
%     A   - constant matrix
%     TOL - tolerance for zeroing (default 1.e-10)
%
%   Outputs: 
%     A0 - transformed matrix
%
%   See also STRIPLZ.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check input data
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1.e-10; end;       
        if ~isnumeric(A), 
           error('Function "findzero" can not be used for this object'); 
        end;
%------------------------------------------------------
%       Zeroing
%------------------------------------------------------
        if norm(A) > tol 
          ind = find(abs(A) < tol*norm(A));
          A(ind) = zeros(size(ind));
        else
          A = A*0; 
        end;
        
%------- End of FINDZERO.M --------- KYuP ----------           
