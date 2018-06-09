function [K,err] = ch2 ( sys, method )
%CH2 H2-optimal continuous-time controller.
%
%     [K,ERR] = CH2 ( SYS, METHOD )
%
%   Inputs:
%     SYS - an LTI system
%     METHOD - calculation method of:
%       'pol' - polynomial (default)
%        'ss' - state-space
%        'wh' - classical Wiener-Hopf
%
%   Outputs:
%     K   - H2-optimal controller 
%     ERR - cost function for the optimal system
%
%   See also H2REG, HINFREG, SDH2, SDL2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        [nout,nin] = size(sys);
        o2 = 1; i2 = 1;
        i1 = nin - i2;
        o1 = nout - o2;
        if i1 < 1, error('Incorrect number of inputs'); end;
        if o1 < 1, error('Incorrect number of outputs'); end;
        if ~exist('method','var'), method = 'pol'; end; 
%------------------------------------------------------
%       State-space method
%------------------------------------------------------
        if isequal(method,'ss')
           sys = ss(sys); 
           K = minreal( zpk(h2reg(sys, 1, 1)) );
           Fcl = minreal(lft(sys,K), sqrt(eps),0);
           err = norm(Fcl)^2;
           err = sqrt(err);
           return;
        end;
%------------------------------------------------------
%       Transform to negative feedback
%------------------------------------------------------
        sys = zpk(sys);
        sys(end,:) = - sys(end,:);
%------------------------------------------------------
%       Extract 4 blocks
%------------------------------------------------------
        P11 = sys(1:o1,1:i1);
        P12 = sys(1:o1,i1+1:nin);
        P21 = sys(o1+1:nout,1:i1);
        P22 = sys(o1+1:nout,i1+1:nin);
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        A0 = minreal(P12'*P12);
        A1 = minreal(P21*P21');
        A = A0 * A1;
        B = minreal(P21*P11'*P12);
        %Edited by RCG to avoid the use of trace
        PT = (P11'*P11);
        n = size(PT,1);
        E = PT(1,1);
        for jj=2:n
            E = E + PT(jj,jj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        E = minreal(E); 
%------------------------------------------------------
%       Solve general quadratic problem
%------------------------------------------------------
        PCancel = intpoles ( sys );  
        err = 0;
        if isequal(method,'pol')
          if nargout < 2 
               K = polquad ( A, B, E, P22, PCancel );
          else [K,err] = polquad ( A, B, E, P22, PCancel );
          end;
        elseif isequal(method,'wh')
           [K,err] = whquad  ( A, B, E, P22, PCancel ); 
        end;
        err = sqrt(err);
        
%------- End of CH2.M --------- KYuP ----------           
        
