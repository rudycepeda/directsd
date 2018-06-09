function F = z2zeta ( F, tol )
%Z2ZETA Change variable z -> 1/z in tf-model.
%       
%     FZ = Z2ZETA ( F, TOL ) 
%
%   Inputs:
%     F(z) - discrete LTI model as tranmsfer function 
%     TOL  - tolerance for determination of zero roots 
%            (default sqrt(eps)) 
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
        [N,D] = tfdata ( F );      
        for r=1:rows
           for c=1:cols
               n = striplz(N{r,c}, tol);
               d = striplz(D{r,c}, tol);
               N{r,c} = striplz(fliplr(n), tol);
               D{r,c} = striplz(fliplr(d), tol);
           end;
        end;        
%--------------------------------------------------
%           Final result
%--------------------------------------------------
        set ( F, 'num', N, 'den', D );

%------- End of Z2ZETA.M --------- KYuP ----------           
