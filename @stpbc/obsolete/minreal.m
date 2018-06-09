function sys = minreal ( sys, recFlag )
%MINREAL Minimal realization for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine system sizes
%------------------------------------------------------
    n = size(sys.a,1);
    m = size(sys.b,2);    
    p = size(sys.c,1);    
    if ~exist('recFlag','var'), recFlag = 0; end;
%------------------------------------------------------
%       Eliminate explicit uncontrollable modes
%------------------------------------------------------
    k = n;
   %------------------------------------------------------
   %  Find the number of zero last rows of B
   %------------------------------------------------------
    while k > 0
      if isequal(sys.b(n-k+1:n,:),sym(zeros(k,m))), break; end;
      k = k - 1;
    end;
    if k == n, % completely uncontrollable system
       sys.a = []; sys.om = []; sys.ups = [];
       sys.b = zeros(0,m); sys.c = zeros(p,0);  
    end;
   %------------------------------------------------------
   %  Explore the structure of A
   %------------------------------------------------------
    while k > 0
      rind = n-k+1:n;
      cind = 1:n-k;
      Ok = zeros(k,n-k);  
      if isequal(sys.a(rind,cind),sym(Ok))
       %-------------------------------------------------
       %  Find suitable transformation for Om and Ups
       %-------------------------------------------------
        if ~(isequal(sys.om(rind,cind),Ok) && ...
             isequal(sys.ups(rind,cind),Ok))
          R = [sys.om(:,cind)'; sys.ups(:,cind)'];  
          N = null(R,'r');
          if size(N,2) < k, continue; end;
          L = null(N','r');
          P = [L'; N'];
          sys.om = P*sys.om;
          sys.ups = P*sys.ups;
        end;  
       %-------------------------------------------------
       %  Perform reduction
       %-------------------------------------------------
        rest = cind;   
        sys.a = sys.a(rest,rest);
        sys.om = sys.om(rest,rest);
        sys.ups = sys.ups(rest,rest);
        sys.b = sys.b(rest,:);
        sys.c = sys.c(:,rest);
        n = n - k;
        break;
      end;
      k = k - 1;
    end;
%------------------------------------------------------
%       Eliminate explicit unobservable modes
%------------------------------------------------------
    if ~recFlag
      sys = minreal(sys.', 1).';
    end;
    
%------- End of MINREAL.M --------- KYuP ----------
