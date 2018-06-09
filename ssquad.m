function [K,err] = ssquad ( A0, A1, B, E, D22 )
%SSQUAD State-space minimization of frequency-domain quadratic functionals.
%
%     [K,ERR] = SSQUAD ( A0, A1, B, E, D22 )
%
%   Inputs:
%     A0(d), A1(d), B(d), E(d) - discrete-time functions
%     D22(d) - discrete model of plant in the loop
%
%   Outputs:
%     K   - optimal discrete controller
%     ERR - cost for the optimal system 
%
%   See also POLQUAD, WHQUAD, QUADERR, SDL2, SDH2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Full standard system
%------------------------------------------------------
        [dsys,X] = abe2std ( A0, A1, B, E, D22 );
        K = h2reg(dsys, 1, 1);
%------------------------------------------------------
%       Cost value (without truncated parts!)
%------------------------------------------------------
        if nargout > 1
          Fcl = minreal(lft(dsys, K), sqrt(eps), 0);
          err = norm(Fcl)^2 + intz2b(X);
        end;
        K = minreal(zpk(K), 1e-3);

%------- End of SSQUAD.M --------- KYuP ----------           
        
