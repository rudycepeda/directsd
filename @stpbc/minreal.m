function [sys,Q,P,k] = minreal ( sys0, recFlag )
%MINREAL Minimal realization for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
    if ~exist('recFlag','var'), recFlag = 0; end;
    sys = sys0;    
    isSym = isa(sys.a,'sym');
    n = size(sys.a,1);
    Q = eye(n); P = Q;
%------------------------------------------------------
%       Construct controllability matrix
%------------------------------------------------------
    Co = ctrb(sys);
    N = null(Co.').';
%------------------------------------------------------
%  Remove uncontrollable modes
%------------------------------------------------------
    kCo = size(N,1);
    while kCo > 0 
      if isSym
           N = double(N);  
           L = null(N,'r').';
      else L = null(N).'; end;
      Qc = [L; N]; Pc = Qc;
      sys = eqtrans ( sys, Qc, Pc );
     %-------------------------------------------------
     %  Block (rind,cind) of A, Om and Ups should be zero
     %-------------------------------------------------
      rind = n-kCo+1:n;
      cind = 1:n-kCo;
     %-------------------------------------------------
     %  Find suitable transformation for Om and Ups
     %-------------------------------------------------
      if isSym
        Ok = zeros(kCo,n-kCo);  
        OmUpsOK = isequal(sys.om(rind,cind),Ok) && ...
                  isequal(sys.ups(rind,cind),Ok);  
      else
        OmUpsOK = norm(sys.om(rind,cind)) < 1e-3 && ...
                  norm(sys.ups(rind,cind)) < 1e-3;  
      end;
      Pa = eye(n);      
      if ~OmUpsOK
        R = double([sys.om(:,cind)'; sys.ups(:,cind)']);  
        N = null(R,'r');
        if size(N,2) >= kCo, 
          L = null(N','r');
          Pa = [L'; N'];
          sys.om = Pa*sys.om;
          sys.ups = Pa*sys.ups;
          OmUpsOK = 1;
        end;
      end;  
     %----------------------------------------------
     %  Perform reduction
     %----------------------------------------------
      if OmUpsOK  
        rest = cind;   
        sys.a = sys.a(rest,rest);
        sys.om = sys.om(rest,rest);
        sys.ups = sys.ups(rest,rest);
        sys.b = sys.b(rest,:);
        sys.c = sys.c(:,rest);
        Q = Qc;
        P = Pa*Pc;
        n = n - kCo;
        break;
      else
        sys = sys0;
        if isSym  
          kCo = 0;  
        else
          kCo = kCo - 1;  
          N = nullr(Co.',kCo).';          
        end;
      end;
    end;    
%------------------------------------------------------
%       Eliminate unobservable modes
%------------------------------------------------------
    if ~recFlag
      [sysOb,Qo,Po,kOb] = minreal(sys.', 1);
      perm = [n-kOb+1:n 1:n-kOb];
      Qo = Qo(perm,:); Po = Po(perm,:);
      sys = sysOb.';
      Qo = blckdiag(inv(Qo.'),eye(kCo)); 
      Po = blckdiag(inv(Po.'),eye(kCo));
      Q = Qo*Q;
      P = Po*P;
      k = kCo + kOb;
    else
      k = kCo;  
    end;
    
%------- End of MINREAL.M --------- KYuP ----------

function N = nullr ( A, k )
    [u,s,v] = svd ( A );
    c = size(A,2); c0 = c - k;
    ind = c0+1:c;
    s(:,ind) = 0;
    N = v*[zeros(c0,k); eye(k)];

%------- End of NULL.M --------- KYuP ----------
