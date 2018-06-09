function C = ctrb ( sys )
%CTRB Controllability matrix for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine system sizes
%------------------------------------------------------
    n = size(sys.a,1);
    m = size(sys.b,2);
%------------------------------------------------------
%       Contsruct controllability matrix
%------------------------------------------------------
    Ai = sys.a;
    C = [sys.b zeros(n,m*(n-1))];
    for i=1:n-1
      ind = (1:m)+m*i;
      C(:,ind) = Ai*sys.b;
      Ai = Ai*sys.a;
    end;
    
%------- End of CTRB.M --------- KYuP ----------
