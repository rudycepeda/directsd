function F = z2zeta ( F, tol )
%Z2ZETA Change variable z -> 1/z in zpk-model.
%       
%     FZ = Z2ZETA ( F, TOL ) 
%
%   Inputs:
%     F(z) - zero-pole-gain form of a discrete LTI model
%     TOL  - tolerance for determination of zero roots (default sqrt(eps)) 
%
%   Outputs:
%     FZ(z) - LTI-model with inverted variable
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~isdt(F), error('Change 1/z for z is applicable for discrete models only'); end;
        if ~exist('tol', 'var'), tol = sqrt(eps); end;       
%------------------------------------------------------
%       Loop over all items
%------------------------------------------------------
        [rows,cols] = size(F);
        [Z,P,K] = zpkdata ( F );      
        for r=1:rows
           for c=1:cols
               zi = Z{r,c};
               pi = P{r,c};               
             %--------------------------------------------
             %   Separate poles and zeros at 0
             %--------------------------------------------
               ind = find(abs(zi) < tol);
               nz0 = length(ind);
               zi(ind) = [];
               ind = find(abs(pi) < tol);
               np0 = length(ind);
               pi(ind) = [];               
             %--------------------------------------------
             %   Inverse remaining poles and zeros               
             %--------------------------------------------
               K(r,c) = real ( K(r,c) * prod(-zi) / prod(-pi) );
               if ~isempty(zi), zi = 1 ./ zi; end;
               if ~isempty(pi), pi = 1 ./ pi; end;
               n0 = length(pi) - length(zi);
               np0 = np0 + n0;
               if nz0 >= np0, pi = [pi; zeros(nz0-np0,1)];
               else           zi = [zi; zeros(np0-nz0,1)]; end;                
             %--------------------------------------------
             %   Return data back
             %--------------------------------------------              
               Z{r,c} = zi;
               P{r,c} = pi;               
           end;
        end;        
%--------------------------------------------------
%           Final result
%--------------------------------------------------
        set ( F, 'z', Z, 'p', P, 'k', K );

%------- End of Z2ZETA.M --------- KYuP ----------           
