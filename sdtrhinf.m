function [K,err] = sdtrhinf ( sys, T, H )
%SDTRHINF AHinf-optimal controller for sampled-data tracking systems.
%
%     [K,ERR] = SDTRHINF ( SYS, T, H )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%     K   - AHinf-optimal discrete controller
%     ERR - cost function for the optimal system
%
%   See also SDL2, SDTRHINFERR, SDAHINF, POLHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        if ~exist('H', 'var'), H = []; end;        
        sys = zpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Construct D22
%------------------------------------------------------
        P22 = sys(end,end);
        D22 = dtfm(P22, T, 0, H)';
        D22 = setpoles(D22, exp(-allpoles(P22)*T) );
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A,B,E] = sdl2coef ( sys, T, H );
%------------------------------------------------------
%       Solve Hinf-type problem
%------------------------------------------------------
        [K,err] = polhinf ( A, B, E, D22 );
        err = sqrt(T)*err;
        
%------- End of SDTRHINF.M --------- KYuP ----------           
        
