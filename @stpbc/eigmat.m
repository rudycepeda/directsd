function [M,eigA,eAT] = eigmat ( sys )
%EIGMAT Test matrix, which should be singular for eigenvalues.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        eigA = sys.a + sys.b*sys.c/sym('gamma^2');
        sys.a = eigA;    
        [M,eAT] = xi ( sys ); 

%------- End of EIGMAT.M --------- KYuP ----------
