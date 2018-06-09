function C = obsv ( sys )
%OBSV Observability matrix for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine system sizes
%------------------------------------------------------
    n = size(sys.a,1);
    p = size(sys.c,1);
%------------------------------------------------------
%       Contsruct observability matrix
%------------------------------------------------------
    Ai = sys.a;
    C = [sys.c; zeros(p*(n-1),n)];
    for i=1:n-1
      ind = (1:p)+p*i;
      C(ind,:) = sys.c*Ai;
      Ai = Ai*sys.a;
    end;
    
%------- End of OBSV.M --------- KYuP ----------
