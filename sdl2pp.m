function [K,err] = sdl2 ( sys, T, H, poles )
%SDL2PP L2-optimal pole-placement controller for sampled-data systems.
%
%     [K,ERR] = SDL2PP ( SYS, T, P )
%     [K,ERR] = SDL2PP ( SYS, T, H, P )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     H   - generalized hold as an LTI model 
%     P   - desired closed-loop poles in 'z'
%
%   Outputs:
%     K    - L2-optimal pole-placement controller as an LTI system
%     ERR  - cost function for the optimal system
%
%   See also SDL2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 26-Dec-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('H','var'), H = []; 
        elseif ~isa(H,'lti')
          poles = H; H = [];  
        end;
        poles = poles(:);
        sys = zpk(sys);
%------------------------------------------------------
%       Obtain characteristic polynomial in 'zeta'
%------------------------------------------------------
        poles = poles(find(abs(poles) > 1e-3));
        Delta = poly(1./poles);
%------------------------------------------------------
%       Check static gain
%------------------------------------------------------
        data.dK0 = 1; % no additional terms
        [gKsi,data.gK,aDelta,bDelta] = psigain ( sys, T, Delta );
        pStatic = ~isnan(gKsi);
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Construct D22
%------------------------------------------------------
        P22 = minreal(sys(end,end));
        poles22 = allpoles(P22);
        plant = dtfm(P22, T, 0, H);
        plant = setpoles(plant, exp(poles22*T) );
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A,B] = sdl2coef ( sys, T );
%------------------------------------------------------
%       Find controller order
%------------------------------------------------------
        p = deg(P22);
        nAll = length(poles);
        ell = nAll - p;
        if ell < p-1+pStatic
          warning('Too few poles, %d poles at zero added',p-1-ell);
          poles = [zeros(p-1-ell,1); poles];
          ell = p-1+pStatic;
        end;
        degKsi = ell - p - pStatic;
%------------------------------------------------------
%       Constrict data structure
%------------------------------------------------------
        data.A = A; data.B = B; 
        data.plant = plant;
%------------------------------------------------------
%       Solve optimization problem
%------------------------------------------------------
        K = optkdelta ( data, Delta, degKsi, gKsi );
%------------------------------------------------------
%       Find cost function
%------------------------------------------------------
        err = sdl2err ( sys, -K );

%------- End of SDL2PP.M --------- KYuP ----------           
        
