function sys = inv ( sys )
%INV Inversion for STPBC B = A^{-1}.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    if isnumeric(double(sys.d))
         Dinv = inv(sys.d);
    else Dinv = sym(['inv(' char(sys.d) ')']); end;
    sys.a = sys.a - sys.b*Dinv*sys.c;
    sys.b = sys.b*Dinv;
    sys.c = -Dinv*sys.c;
    sys.d = Dinv;

%------- End of MTIMES.M --------- KYuP ----------
