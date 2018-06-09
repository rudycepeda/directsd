function err = sd2doferr ( sys, K, KR, H )
%SD2DOFERR Integral quadratic error in 2-DOF sampled-data system.
%
%     ERR = SD2DOFERR ( SYS, K, KR, H )
%
%   Inputs:
%     SYS - an LTI system
%     K   - discrete-time loop controller
%     KR  - discrete-time reference controller
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%     ERR - L2-norm of continuous-time error
%
%   See also SDL2ERR, SD2DOF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('H','var'), H = []; end;
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys = zpk(sys);
        sys(end-1,:) = -sys(end-1,:);
        K = - K;
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        T = K.Ts;
        [A,B,E,kCoef,n,d] = sd2dofcoef ( sys, T, H );
%------------------------------------------------------
%       Construct b and Delta
%------------------------------------------------------
        K = K'; KR = KR';
        [a,b] = tf2nd ( K );        
        Delta = a*n + b*d;
%------------------------------------------------------
%       Construct PsiX = (k*a + KR*b)/Delta
%------------------------------------------------------
        PsiX = (zpk(kCoef*a,1,T) + KR*zpk(b,1,T))*zpk(1,Delta,T);  
        PsiX = minreal(PsiX);
%------------------------------------------------------
%       Compute error
%------------------------------------------------------
        err = quaderr ( A, B, E, PsiX );
        
%------- End of SD2DOFERR.M --------- KYuP ----------           
        
