function err = sdahinorm ( sys, K, H )
%SDAHINORM Associated Hinf-norm for sampled-data system.
%
%     ERR = SDAHINORM ( SYS, K, H )
%
%   Inputs:
%     SYS - an LTI system
%     K   - discrete-time controller
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%     ERR - AHinf-norm of sampled-data system
%
%   See also SDAHINF, SDHINORM, SDTRHINFERR, SDH2NORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
        K = - K;
        if ~exist('H','var') || isempty(H)
             err = sdhinferr ( sys, K, 'sdh2coef' );
        else err = sdhinferr ( sys, K, 'sdh2coef', H ); 
        end;
%------------------------------------------------------
%       Special scaling for compatibility with Hinf-norm
%------------------------------------------------------
        err = sqrt(K.Ts)*err;

%------- End of SDAHINORM.M --------- KYuP ----------           
        
