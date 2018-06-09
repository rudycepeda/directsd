function sys = uminus ( sys )
%UMINUS Unary minus for STPBC B = - A.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        sys.c = - sys.c;
        sys.d = - sys.d;

%------- End of UMINUS.M --------- KYuP ----------
