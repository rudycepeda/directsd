function [V,err] = intz ( varargin )
%INTZ Compute contour integral using direct calculation.
% 
%     [V,ERR] = INTZ ( X1, X2, ...  )
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
%       Main loop
%------------------------------------------------------
        V = 0; P = 0;
        for i=1:nargin
          S = zpk(varargin{i});
         %--------------------------------------------- 
         %   Remove poles at z=1
         %--------------------------------------------- 
          pp = S.p{1};
          ind = find(abs(abs(pp)-1) < 1e-4);
          if ~isempty(ind),
            [P0,S] = separss(S, pp(ind));  
            P = sumzpk(P, P0);
          end;
         %--------------------------------------------- 
         %   Make the function positive definite
         %--------------------------------------------- 
          a = 2*norm(S,'inf');
          F = sfactor(S + a);
          V = V + norm(F')^2 - a;
        end;
%------------------------------------------------------
%       Truncation error
%------------------------------------------------------
        err = norm(tfdata(P,'v'));

%------- End of INTZ2B.M --------- KYuP ----------           
