function err = sderr ( sys, K, fname, t, H )
%SDERR Squared integral error for sampled-data system.
%
%     ERR = SDERR ( SYS, K, FNAME, T )
%     ERR = SDERR ( SYS, K, FNAME, H )
%     ERR = SDERR ( SYS, K, FNAME, T, H )
%
%   Inputs:
%     SYS   - an LTI system
%     K     - discrete-time controller
%     FNAME - name of function for computing coefficients
%     T     - time instant for special H2-problems
%     H     - generalized hold as an LTI model 
%
%   Outputs:
%     ERR   - L2-norm of continuous-time error
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('t','var'), t = []; 
        elseif isa(t,'lti')
           if exist('H','var'),
              error('Parameter T must be a constant'); 
           end;
           H = t; t = [];
        end;
        if ~exist('H', 'var'), 
           H = ss(0, 1, 1, 0); % scalar ZOH by default
        end;
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        sys = zpk(sys);
        P22 = sys(end,end);
        T = K.Ts;
        switch fname
          case 'sdh2coef'  
            if isempty(t),
                 cmd = sprintf('[A,B,E] = %s ( sys, T, [], H );', fname);             
            else cmd = sprintf('[A,B,E] = %s ( sys, T, t, H );', fname); 
            end;        
          case 'sdl2coef'  
            if ~isempty(t), error('Parameter ''t'' must be empty'); end;
            cmd = sprintf('[A,B,E] = %s ( sys, T, H );', fname);             
          otherwise
            error('Unknown function name ''%s''', fname); 
        end;
        eval ( cmd );
%------------------------------------------------------
%       Construct D22 and M
%------------------------------------------------------
        D22 = dtfm(P22, T, 0, H);
        M = zpk(K / (1 + K*D22)); %M = zpk(feedback ( K, D22 ));
        ind = find(abs(M.z{1}) < sqrt(eps));
        M.z{1}(ind) = 0;
        ind = find(abs(M.p{1}) < sqrt(eps));
        M.p{1}(ind) = 0;
        M = M';
%------------------------------------------------------
%       Compute error
%------------------------------------------------------
        err = quaderr ( A, B, E, M );            %sqrt(err)
        %err = quaderrss ( A0, A1, B, E, D22', K );% sqrt(err)
        %err = quaderr2  ( A0, A1, B, E, D22, K ); % sqrt(err)

%------- End of SDERR.M --------- KYuP ----------           
        
