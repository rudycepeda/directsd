function [V,err] = intz2b ( varargin )
%INTZ2B Compute contour integral using bilinear transform.
% 
%     [V,ERR] = INTZ2B ( X1, X2, ...  )
%
%   Inputs:
%     X1, X2, ... - symmetric discrete-time functions
%
%   Outputs:
%     V = int_Gamma (X1 + X2 + ...) dz/z
%     ERR - truncation error (norm of the polynomial caused
%           by poles at z=1)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Bilinear transform parameters
%------------------------------------------------------
        b = 1;
        BL = fbilin ( b );
%------------------------------------------------------
%       Main loop
%------------------------------------------------------
        V = 0; P = 0;
        for i=1:nargin
          S = sdzpk(varargin{i});          
         %--------------------------------------------- 
         %   Make the function positive definite
         %--------------------------------------------- 
          Sb  = bilintr(S, b);
          [P0,S0] = improper ( Sb, 'symm' );
         %--------------------------------------------- 
         %   Correct improper part
         %--------------------------------------------- 
          if norm(P0) > sqrt(eps)  
            val1 = polyval(P0, 1);
            P0 = P0 - val1;
            S0 = S0 + val1;
            P = P0+P;
          end;
         %--------------------------------------------- 
         %   Make the function positive definite
         %--------------------------------------------- 
          a = 2*norm(S0,'inf');
          S0 = (S0 + a)*BL;
          F = sfactor(S0);
          V = V + norm(F)^2 - a;
        end;
%------------------------------------------------------
%       Truncation error
%------------------------------------------------------
        err = norm(P);

%------- End of INTZ2B.M --------- KYuP ----------           
