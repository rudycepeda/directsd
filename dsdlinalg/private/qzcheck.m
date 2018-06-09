function flag = qzcheck ( A0, A, Q, Z, tol )
%QZCHECK Check correctness of QZ-transformation.
%         
%     QZCHECK ( A0, A, Q, Z )
%     QZCHECK ( A0, A, Q, Z, TOL )
%     F = QZCHECK ( A0, A, Q, Z )
%     F = QZCHECK ( A0, A, Q, Z, TOL )
%
%   Inputs:
%     A0, A - matrices of the same size
%     Q, Z  - unitary matrices such that A = Q*A0*Z
%
%   Outputs: 
%     F = 1, if norm(Q*A0*Z-A) < TOL*norm(A) 
%     F = 0, otherwise
%
%   If no output parameters are given, reports error on F = 0.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1.e-6; end;       
        flag = norm(Q*A0*Z - A) < tol*norm(A);
        if ~flag  &&  nargout == 0
           error('Incorrect QZ-transformation');     
        end;

%------- End of QZCHECK.M --------- KYuP ----------           
