function [K,err,errPrev] = sdl2 ( sys, T, H, method )
%SDL2 L2-optimal controller for sampled-data systems.
%
%     [K,ERR,PERR] = SDL2 ( SYS, T, H, METHOD )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     H   - generalized hold as an LTI model 
%     METHOD - calculation method of:
%        'pol' polynomial (default)
%        'ss'  state-space (lifting)
%        'ssf' state-space based on frequency domain
%        'wh'  classical Wiener-Hopf
%
%   Outputs:
%     K    - L2-optimal discrete controller as an LTI system
%     ERR  - cost function for the optimal system
%     PERR - limiting cost function for the optimal systems
%            as preview interval tends to Inf
%
%   See also SDH2, SDH2HINF.

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
        if ~exist('method','var'), method = 'ss'; end;        
%------------------------------------------------------
%       Check method
%------------------------------------------------------
        td = totaldelay(sys);
        noDelay = all(td(:) < eps);
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
        if ~isequal(method,'pol') &&  ~isequal(method,'ss') && ...
           ~isequal(method,'ssf') && ~isequal(method,'wh')
          error('Unknown calculation method ''%s''',method);  
        end;
%------------------------------------------------------
%       State-space method
%------------------------------------------------------
        if isequal(method,'ss')
           dsysL2 = sdh2simple ( sys, T );
           K = sdzpk(h2reg ( dsysL2, 1, 1 ));
           K = minreal(K,1e-3);
           Fcl = minreal(lft(dsysL2,K), sqrt(eps),0);
           err = norm(Fcl)^2;
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
        D22 = dtfm(P22, T, 0, H)';
        D22 = setpoles(D22, exp(-allpoles(P22)*T) );
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        [A,B,E,A0,A1] = sdl2coef ( sys, T );
%------------------------------------------------------
%       Solve general quadratic problem
%------------------------------------------------------
        if isequal(method,'pol')
           switch nargout 
            case {0, 1}, 
              K = polquad ( A, B, E, D22 );
            case 2,      
              [K,err] = polquad ( A, B, E, D22 );
            case 3,
              [K,err,err0,err2] = polquad ( A, B, E, D22 );
              errPrev = err0 + err2;   
           end;
        elseif isequal(method,'wh')
           if nargout < 3, 
             [K,err] = whquad  ( A, B, E, D22 ); 
           else
             [K,err,err0,err2] = whquad  ( A, B, E, D22 ); 
             errPrev = err0 + err2;   
           end;
        elseif isequal(method,'ssf')
           errPrev = 0;  
           [K,err] = ssquad ( A0, A1, B, E, D22 );
        end;

%------- End of SDL2.M --------- KYuP ----------           
        
