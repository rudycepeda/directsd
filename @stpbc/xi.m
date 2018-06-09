function [Xi,eAT,invXi] = xi ( sys, T )
%XI Xi matrix for STPBC.
%

%------------------------------------------------------
% Copyright 1999-2007 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine system sizes
%------------------------------------------------------
    n = size(sys.a,1);
    isSym = isa(sys.a,'sym');
%------------------------------------------------------
%       Build exp(AT) using structure
%------------------------------------------------------
    if isSym
      Ar = rand(n) .* sign(sys.a);
      mask = sign(expm(Ar));
      eAT = sym(zeros(n));
      for i=1:n, 
        for j=1:n
          if mask(i,j),
             eAT(i,j) = sym(sprintf('G%d%d',i,j)); 
          end;
        end;
      end;
%------------------------------------------------------
%       Build exp(AT) using numerical exponential
%------------------------------------------------------
    else
      if ~exist('T','var'), T = 1; end;
      eAT = expm(sys.a*T);
    end;
%------------------------------------------------------
%       Xi matrix
%------------------------------------------------------
    Xi = sys.om + sys.ups*eAT;
%------------------------------------------------------
%       Inverted Xi matrix
%------------------------------------------------------
    if nargout > 2
      if isSym
           invXi = invmat(Xi);
      else invXi = inv(Xi); end;        
    end;

%------- End of XI.M --------- KYuP ----------

%######################################################
function sgn = sign ( x, tol )
    if ~exist('tol','var'), tol = 1e-6; end;
    [rr,cc] = size(x);
    sgn = ones(rr,cc);  
    if isa(x,'sym'), 
         ind0 = find(x == 0); 
    else ind0 = find(abs(x) < tol*norm(x)); end;
    sgn(ind0) = 0;

%------- End of SIGN.M --------- KYuP ----------

    