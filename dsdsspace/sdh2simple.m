function [dsys,DP11,DP12,DP21,DP22] = sdh2simple ( csys, T, o2, i2, tau )
%SDH2SIMPLE Simple H2-equivalent model of a sampled-data system.
%
%     [DSYS,DP11,DP12,DP21,DP22] = SDH2SIMPLE ( SYS, T, NMEAS, NCON, TAU )
%
%   Inputs:
%     SYS - continuous-time standard system 
%               | z | = | P11 P12 | | w |
%               | y |   | P21 P22 | | u |
%  	  T   - sampling period
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1)
%     TAU   - input time-delay (default 0)
%
%   Outputs:
%	  DSYS  - discrete-time H2-equivalent of the system
%               | zd | = | DP11 DP12 | | dw |
%               | yd |   | DP21 DP22 | | du |
%
%   See also SDFAST, SDGH2MOD, SDH2.

% References:
%   [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%       Springer-Verlag, Berlin etc., 1995.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 29-Oct-2006 $
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
	    sys = minreal(ss(csys), sqrt(eps), 0);
	    if ~isct(sys)
	       error('Initial system must be a continuous time model');
	    end;
        if ~exist('tau','var'), tau = 0; end;
        if tau > T, warning('Tau > T, assuming tau=0.'); end;
        [a,b,c,d] = ssdata(sys);
        n = size(a,1);
        [nout,nin] = size(sys);        
%------------------------------------------------------
%	Check validity
%------------------------------------------------------
        if ~exist('i2','var'), i2 = 1; end;
        if ~exist('o2','var'), o2 = 1; end;
	    i1 = nin - i2;
	    o1 = nout - o2;
	    if i1 < 1, error('Negative or zero number of reference inputs'); end;
	    if o1 < 1, error('Negative or zero number of working outputs');	end;
	    if i1 >= nin, error('No control signals specified'); end;
	    if o1 >= nout, error('No measured signals specified'); end;
%------------------------------------------------------
%	Extract partial matrices
%	ans check assumptions
%------------------------------------------------------
	    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
	    if norm(d11) > eps, warning('D11 must be (assumed) zero');	end;
	    if norm(d21) > eps, warning('D21 must be (assumed) zero');	end;
	    if norm(d22) > eps, warning('D22 must be (assumed) zero');	end;
	    %d11 = 0*d11;
	    %d21 = 0*d21;
	    %d22 = 0*d22;
%------------------------------------------------------
%	Loop discretization
%------------------------------------------------------
        Abar = [   a   b2
                zeros(i2,n+i2)];    
        E = expm(Abar*T);
        Ad  = E(1:n,1:n);
        B2d = E(1:n,n+1:end);        
%------------------------------------------------------
%	Output discretization
%------------------------------------------------------
        c1d12 = [c1 d12];
        Q = c1d12'*c1d12;
        Qd = intaxa(Abar', Q, T);
        C1dD12d = rectchol(Qd);
        C1d  = C1dD12d(:,1:n);
        D12d = C1dD12d(:,n+1:end);
%------------------------------------------------------
%	Final discrete model with the same H2-optimal controller
%------------------------------------------------------   
        i1d = i1;
        o1d = size(C1d,1);
        D11d = zeros(o1d,i1d);
        if tau == 0 
           B1d = b1;
        else
           B1d = expm(a*(T-tau))*b1;
           Z = intaxa(a, b1*b1', T-tau);  
           gamma = trace(c1*Z*c1')
           D11d(1,1) = sqrt(gamma);
        end;
        Bd = [B1d B2d];
        Cd = [C1d; c2];
        Dd = [   D11d   D12d  
               zeros(o2,i1d+i2) ];
        dsys = minreal ( ss ( Ad, Bd, Cd, Dd, T ), sqrt(eps), 0 );       
%------------------------------------------------------
%	Additional data if necessary
%------------------------------------------------------   
        if nargout > 1
           DP11 = minreal( ss(Ad, B1d, C1d, zeros(o1d,i1d), T), sqrt(eps), 0 );
        end;    
        if nargout > 2
           DP12 = minreal( ss(Ad, B2d, C1d, D12d, T), sqrt(eps), 0 );
        end;
        if nargout > 3
           DP21 = minreal( ss(Ad, B1d, c2, zeros(o2,i1d), T), sqrt(eps), 0 );
        end;
        if nargout > 4
           DP22 = minreal( ss(Ad, B2d, c2, zeros(o2,i2), T), sqrt(eps), 0 );            
        end;      
    
%------- End of SDH2SIMPLE.M --------- KYuP ----------
