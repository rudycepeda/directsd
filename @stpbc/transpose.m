function sys = transpose(sys)
%TRANSPOSE Transpose B = A.'
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    n = size(sys.a,1);
    sys.a = transpsym(sys.a);
    Btemp = sys.b;
    sys.b = transpsym(sys.c);
    sys.c = transpsym(Btemp);
    sys.d = transpsym(sys.d);
    sys.om = sys.om';
    sys.ups = sys.ups';

%------- End of TRANSPOSE.M --------- KYuP ----------

