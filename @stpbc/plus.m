function sys = plus ( sys1, sys2 )
%PLUS Sum of two STPBC C = A + B.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    n1 = size(sys1.a,1);
    n2 = size(sys2.a,1);
    sys = stpbc();
    sys.a = blckdiag ( sys1.a, sys2.a );
    sys.om = blckdiag ( sys1.om, sys2.om );
    sys.ups = blckdiag ( sys1.ups, sys2.ups );
    sys.b = [sys1.b; sys2.b];
    sys.c = [sys1.c sys2.c];
    sys.d = sys1.d + sys2.d;

%------- End of PLUS.M --------- KYuP ----------
