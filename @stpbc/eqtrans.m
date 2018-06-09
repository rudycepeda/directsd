function sys = eqtrans (sys, Q, P)
%EQTRANS Equivalent transformation for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    n = size(sys.a,1);
    if ~exist('P','var'), P = Q; end;
    sys.a = Q*sys.a/Q;
    sys.b = Q*sys.b;
    sys.c = sys.c/Q;
    sys.om = P*sys.om/Q;
    sys.ups = P*sys.ups/Q;    

%------- End of CTRANSPOSE.M --------- KYuP ----------
