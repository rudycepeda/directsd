function err = sdh2norm ( sys, K, t, H, method )
%SDH2NORM H2-norm of sampled-data systems.
%
%     VM = SDH2NORM ( SYS, K )
%     VM = SDH2NORM ( SYS, K, H, METHOD )
%     VM = SDH2NORM ( SYS, K, METHOD )
%     V  = SDH2NORM ( SYS, K, t )
%     V  = SDH2NORM ( SYS, K, t, METHOD )
%     V  = SDH2NORM ( SYS, K, t, H, METHOD )
%
%   Inputs:
%     SYS - an LTI system
%     K   - discrete-time controller
%     t   - time instant for calculating variance
%     H   - generalized hold as an LTI continuous-time system
%     METHOD - calculation method of:
%        'ss'  state-space method (inapplicable for time-delay
%              systems and for a specified t)   
%        'pol' polynomial (default for delayed systems)
%
%   Outputs:
%     VM - H2-norm of continuous-time error
%     V  - standard deviation at specified time instants
%
%   See also SDH2, SDAHINORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 22-Sep-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('t','var'), t = []; 
        elseif isa(t,'lti')
           if exist('H','var')
               method = H; 
           end
           H = t; 
           t = [];
        elseif ischar(t)
           method = t; 
           t = []; 
        end
        if ~exist('H','var')
            H = []; 
        elseif ischar(H)
          method = H; 
          H = [];  
        end
        K = zpk(K);
%------------------------------------------------------
%       Check method
%------------------------------------------------------
        td = totaldelay(sys);
        noDelay = all(td(:) < eps);
        if ~exist('method','var'), 
          if noDelay
               method = 'ss'; 
          else method = 'pol'; 
          end; 
        end;
%------------------------------------------------------
        msgdelay = 'System has delays, state-space computations impossible'; 
        if isequal(method,'ss')
          if ~noDelay, warning(msgdelay);   method = 'pol'; end;
        end;
%------------------------------------------------------
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
           err = sderr ( sys, K, 'sdh2coef', t, H );
           err = sqrt(err);
           return;
        end;
%------------------------------------------------------
%       For delay-free systems, use state=space method
%------------------------------------------------------
        if isempty(t)
          dsysH2 = sdgh2mod ( sys, K.Ts );
          dsysH2 = sdzpk(dsysH2);
          Fcl = minreal(lft(dsysH2,K),sqrt(eps));
          err = norm(Fcl);
        else
          err = t; 
          for i=1:length(t)
            dsysH2 = sdgh2mod ( sys, K.Ts, 1, 1, t(i) );
            dsysH2 = sdzpk(dsysH2);
            Fcl = minreal(lft(dsysH2,K), sqrt(eps));
            err(i) = norm(Fcl);
          end;
        end

%------- End of SDH2NORM.M --------- KYuP ----------           
        
