function [K,err,errPrev] = sdh2 ( sys, T, t, H, method )
%SDH2 H2-optimal controller for sampled-data systems.
%
%     [K,ERR,PERR] = SDH2 ( SYS, T, t, H, METHOD )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     t   - time instant for optimization
%     H   - generalized hold as an LTI model 
%     METHOD - calculation method of:
%       'pol' polynomial (default)
%       'ss'  state-space (lifting)
%       'ssf' state-space based on frequency domain
%       'wh'  classical Wiener-Hopf
%
%   Outputs:
%     K   - H2-optimal discrete controller as an LTI system
%     ERR - cost function for the optimal system
%     PERR - limiting cost function for the optimal systems
%           as preview interval tends to Inf
%
%   See also SDL2, CH2, SDH2HINF. 

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 30-Oct-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('t','var'), t = []; 
        elseif isa(t,'lti')
           H = t; t = [];
        elseif ischar(t)
           method = t; t = []; 
        end;
        if ~exist('H','var'), H = []; 
        elseif ischar(H)
          method = H; H = [];  
        end;
%------------------------------------------------------
%       Check method
%------------------------------------------------------
        td = totaldelay(sys);
        noDelay = all(td(:) < eps);
        if ~exist('method','var'), 
          if noDelay, method = 'ss';
          else        method = 'pol'; end;
        end; 
%------------------------------------------------------
        msgdelay = 'System has delays, state-space computations impossible'; 
        msgprev = 'State-space computations impossible for preview analysis';
        msghold = 'Method ''%s'' is not applicable for generalized hold';
        if isequal(method,'ss')
          if ~noDelay, warning(msgdelay);   method = 'pol'; end;
          if nargout > 2, warning(msgprev); method = 'pol'; end;          
          if ~isempty(H), disp(sprintf(msghold, method)); method = 'pol'; end;  
        end;
%------------------------------------------------------
        if ~isequal(method,'pol')  &&  ~isequal(method,'ss') && ...
           ~isequal(method,'ssf') && ~isequal(method,'wh')
          error('Unknown calculation method ''%s''',method);  
        end;
%------------------------------------------------------
%       State-space method
%------------------------------------------------------
        if isequal(method,'ss')
           if isempty(t) 
               [dsysH2,gamma] = sdgh2mod ( sys, T );
           else
               [dsysH2,gamma] = sdgh2mod ( sys, T, 1, 1, t );
           end
           K = minreal(zpk(h2reg ( dsysH2, 1, 1 )));
           Fcl = minreal(lft(dsysH2,K), sqrt(eps),0);
           err = norm(Fcl)^2 + gamma;
           err = sqrt(err);
           return;
        end;
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys = sdzpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Construct D22
%------------------------------------------------------
        P22 = sys(end,end);
        D22 = sdzpk(dtfm(P22, T, 0, H)');
        D22 = setpoles(D22, exp(-allpoles(sdzpk(P22))*T) );
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A,B,E,A0,A1] = sdh2coef ( sys, T, t, H );
%------------------------------------------------------
%       Solve general quadratic problem
%------------------------------------------------------
        PCancel = intpoles ( sys );  
        err = 0;
        if isequal(method,'pol')
          switch nargout 
            case {0, 1}, 
             K = polquad ( A, B, E, D22, PCancel );
            case 2, 
             [K,err] = polquad ( A, B, E, D22, PCancel );
            case 3,            
             [K,err,err0,err2] = polquad ( A, B, E, D22, PCancel );
             errPrev = sqrt(err0 + err2);   
           end;
        elseif isequal(method,'wh')
           if nargout < 3, 
             [K,err] = whquad  ( A, B, E, D22, PCancel ); 
           else
             [K,err,err0,err2] = whquad  ( A, B, E, D22, PCancel ); 
             errPrev = sqrt(err0 + err2);   
           end;
        elseif isequal(method,'ssf')
           if nargout > 1
             errPrev = 0;
             [K,err] = ssquad ( A0, A1, B, E, D22 );
           else
             err = 0;  
             K = ssquad ( A0, A1, B, E, D22 );
           end;
        end;
        err = sqrt(err);
        
%------- End of SDH2.M --------- KYuP ----------           
        
