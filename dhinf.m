function [K,err] = dhinf ( sys )
%DHINF Polynomial Hinf-optimal discrete-time system design.
%
%     [K,ERR] = DHINF ( SYS )
%
%   Inputs:
%     SYS - a discrete-time LTI system
%
%   Outputs:
%     K   - Hinf-optimal discrete controller as an LTI system
%     ERR - cost function for the optimal system
%
%   See also SDAHINF, SDTRHINF, POLHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Transform to BACKWARD shift and negative feedback
%------------------------------------------------------
        sys = z2zeta ( sys );
        sys = zpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------        
        D22 = sys(end,end);
        [A,B,E] = h2coef ( sys );
%------------------------------------------------------
%       Solve Hinf-type problem
%------------------------------------------------------
        [K,err] = polhinf ( A, B, E, D22 );

%------- End of DHINF.M --------- KYuP ----------           
        
