function v = val ( sys, x )
%VAL Value of transfer function at a point (s or z).
%
%     V = VAL ( SYS, P )
%
%   Inputs:
%     SYS - system model
%     P   - value of independent variable
%
%   Outputs:
%     V  - value of transfer function
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       ZPK-form
%------------------------------------------------------
        if isa(sys,'zpk')
          [zz,pp,kk] = zpkdata ( sys, 'v' );
          v = x;
          for i=1:length(x)
            v(i) = kk * prod(x(i) - zz) / prod(x(i) - pp);
          end;
          return;
        end;
%------------------------------------------------------
%       TF-form
%------------------------------------------------------
        if isa(sys,'tf')
          [nn,dd] = tfdata ( sys, 'v' );
          v = x;
          for i=1:length(x)
            v(i) = polyval(nn,x(i)) / polyval(dd, x(i));
          end
          return;
        end;
%------------------------------------------------------
%       SS-form
%------------------------------------------------------
        if isa(sys,'ss')
          [a,b,c,d] = ssdata ( sys );
          n = size(a,1);
          v = x;
          for i=1:length(x)
            v(i) = c*((x(i)*eye(n) - a)\b) + d;
          end;
          return;
        else
          error('Unsupported class ''%s''',class(sys));  
        end;

%------- End of VAL.M --------- KYuP ----------
