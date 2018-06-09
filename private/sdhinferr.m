function err = sdhinferr ( sys, K, fname, H )
%SDHINFERR AHinf-norm for sampled-data system
%
%     ERR = SDHINFERR ( SYS, K, FNAME )
%     ERR = SDHINFERR ( SYS, K, FNAME, H )
%
%   Inputs:
%     SYS   - an LTI system
%     K     - discrete-time controller
%     FNAME - name of function for computing coefficients
%     H     - generalized hold as an LTI model 
%
%   Outputs:
%     ERR   - AHinf-norm of sampled-data system in
%             deterministic or stochastic sense
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
        if ~exist('H', 'var'), H = []; end;
        sys = zpk(sys);
        P22 = sys(end,end);
        T = K.Ts;
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        switch fname
          case 'h2coef'  
            cmd = sprintf('[A,B,E] = %s ( sys );', fname);             
          case 'sdh2coef'  
            cmd = sprintf('[A,B,E] = %s ( sys, T, [], H );', fname);             
          case 'sdl2coef'  
            cmd = sprintf('[A,B,E] = %s ( sys, T, H );', fname);             
          otherwise
            error('Unknown function name ''%s''', fname); 
        end;
        eval ( cmd );
%------------------------------------------------------
%       Construct D22 and M
%------------------------------------------------------
        if sys.Ts ~= K.Ts
             D22 = dtfm(P22, T, 0, H)';
        else D22 = P22; 
        end;
        M = feedback ( K', D22 );
%------------------------------------------------------
%       Compute maximum of the frequency responce
%------------------------------------------------------
        AMM = minreal(A*M*M');
        BM = minreal(B*M);
        X = sumzpk(AMM, E, -BM, -BM');
        ws = warning; warning off;
        try   X = minreals(X, 1e-3);
        catch X = minreal(X, 1e-3); end;
        warning(ws);
%------------------------------------------------------
%       Eliminate poles at z=1
%------------------------------------------------------
        [nX,dX] = tf2nd ( X );
        ind = find(abs(abs(dX.z)-1) < 1e-6);
        if ~isempty(ind)
          v = poln(dX.z(ind), ['r' dX.var]);
          nX = nX / v;
          dX = dX / v;
          X = zpk(nX, dX, X.Ts);
        end;
%------------------------------------------------------
%       Direct frequency response method
%------------------------------------------------------
        %err = sqrt(max(bode(X)));
        %return
        err = sqrt(zinfnorm(X));

%------- End of SDHINFERR.M --------- KYuP ----------           
        
