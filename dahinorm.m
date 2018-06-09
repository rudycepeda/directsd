function err = dahinorm ( sys, K )
%DAHINORM Associated Hinf-norm for discrete-time system.
%
%     ERR = DAHINORM ( SYS, K )
%
%   Inputs:
%     SYS - discrete-time LTI system in FORWARD shift operator
%     K   - discrete-time controller in FORWARD shift operator
%
%   Outputs:
%     ERR - AHinf-norm of discrete-time system
%
%   See also SDAHINORM, SDTRHINFERR,

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
        K = - K;
        sys = z2zeta ( sdzpk(sys) );
        err = sdhinferr ( sys, K, 'h2coef' );

%------- End of DAHINORM.M --------- KYuP ----------           
        
