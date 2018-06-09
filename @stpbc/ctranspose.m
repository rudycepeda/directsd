function sys = ctranspose(sys)
%CTRANSPOSE Complex conjugate transpose B = A'.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    n = size(sys.a,1);
    isSym = isa(sys.a,'sym');
%------------------------------------------------------
%   Boundary conditions
%------------------------------------------------------
    if isSym
         OK = isequal(sys.om,eye(n));
    else OK = (norm(sys.om-eye(n)) < 1e-6); end;
    if OK
      sys.om = zeros(n);
      sys.ups = eye(n);
    else
      %-------------------------------------------------   
      %  Find Om*Ups1 = Ups*Om1
      %-------------------------------------------------   
       N = null([sys.om -sys.ups], 'r');
       if size(N,2) ~= n
         error('Incorrect conjugate')  
       end;
       Ups1 = N(1:n,:)';
       Om1 = N(n+1:end,:)';
       sys.om = Ups1;
       sys.ups = Om1;
    end;
%------------------------------------------------------
%   Other transformations
%------------------------------------------------------
    sys.a = - sys.a';
    Btemp = sys.b;
    %    sys.b = sys.c';
    %    sys.c = - Btemp';
    sys.b = - sys.c';
    sys.c = Btemp';
    sys.d = sys.d';    

%------- End of CTRANSPOSE.M --------- KYuP ----------

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

