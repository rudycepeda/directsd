function err = sdl2err ( sys, K, H, method )
%SDL2ERR Integral quadratic error for sampled-data systems.
%
%     ERR = SDL2ERR ( SYS, K, H, METHOD )
%
%   Inputs:
%     SYS - an LTI system
%     K   - discrete-time controller
%     H   - an extrapolator (default is ZOH)
%     METHOD - calculation method of:
%        'pol' (mandatory for systems with delays)
%        'ss'  (default for delay-free case)
%
%   Outputs:
%     ERR - L2-norm of continuous-time error
%
%   See also SDL2, SD2DOF, SD2DOFERR.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 29-Oct-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('H','var'), H = []; 
        elseif ischar(H)
          method = H; H = [];  
        end;
%------------------------------------------------------
%       Choose method
%------------------------------------------------------
        td = totaldelay(sys);
        noDelay = all(td(:) < eps);
        if ~exist('method','var'), 
          if noDelay  
               method = 'ss'; 
          else method = 'pol'; 
          end; 
        end;
        if ~noDelay  &&  isequal(method,'ss')
          warning('System has delays, state-space computations impossible');  
          method = 'pol';
        end;
        if ~isequal(method,'pol')  &&  ~isequal(method,'ss')
          error('Unknown calculation method ''%s''',method);  
        end;
%------------------------------------------------------
%       Extrapolator
%------------------------------------------------------
        if ~isempty(H), 
          if ~isequal(method,'pol'), 
            disp(sprintf('Method ''%s'' is not applicable for generalized hold', method));  
            method = 'pol';  
          end;
        else
          H = ss(0, 1, 1, 0); % scalar ZOH by default
        end;
%------------------------------------------------------
%       General case: use polynomial method
%------------------------------------------------------
        if isequal(method,'pol')
           K = - K;
           err = sderr ( sys, K, 'sdl2coef', [], H );
           return;
        end;
%------------------------------------------------------
%       Delay-free case: use state-space method
%------------------------------------------------------
        dsysL2 = sdh2simple ( sys, K.Ts );
        dcl = minreal(lft(dsysL2, K), sqrt(eps), 0);
        ws = warning; warning off;
        err = norm(dcl)^2;
        if isinf(err)
           dcl = minreal(zpk(dcl), 1e-3);
           err = norm(dcl)^2;
        end;
        warning(ws);

%------- End of SDL2ERR.M --------- KYuP ----------           
        
