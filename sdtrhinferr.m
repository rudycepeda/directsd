function err = sdtrhinferr ( sys, K, H )
%SDTRHINFERR Associated Hinf-cost for sampled-data tracking systems.
%
%     ERR = SDHINFERR ( SYS, K )
%     ERR = SDHINFERR ( SYS, K, H )
%
%   Inputs:
%     SYS - an LTI system
%     K   - discrete-time controller
%     H   - generalized hold as an LTI model (default ZOH)
%
%   Outputs:
%     ERR - AHinf-norm of sampled-data tracking system
%
%   See also SDTRHINF, SDAHINF, SDAHINORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
        if ~exist('H', 'var'), H = []; end;
        K = - K;
        err = sdhinferr ( sys, K, 'sdl2coef', H );
        err = sqrt(K.Ts)*err;

%------- End of SDHINFERR.M --------- KYuP ----------           
        
