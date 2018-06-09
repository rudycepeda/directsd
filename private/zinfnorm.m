function N = zinfnorm ( X )
%ZINFNORM Compute Linf-norm for discrete-time model.
% 
%     N = ZINFNORM ( X )
%
%   Inputs:
%     X - discrete-time LTI model
%
%   Outputs:
%     N - Linf-norm (maximal value at the unit circle)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Ensure properness
%------------------------------------------------------
        X = minreal(zpk(X));
        nz = length(X.z{1});
        np = length(X.p{1});
        if nz > np, 
          X = X * zpk([], zeros(nz-np,1), 1, X.Ts);  
        end;        
%------------------------------------------------------
%       Computation using state-space algorithm
%------------------------------------------------------
        N = norm(X,'inf');

%------- End of ZINFNORM.M --------- KYuP ----------           
