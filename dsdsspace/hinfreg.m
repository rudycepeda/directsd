function [Kopt,gamma,ClSys] = hinfreg ( sys, o2, i2, options )
%HINFREG Hinf-optimal controller for an LTI-system.
%
%     [K,GAMMA,CL] = HINFREG ( SYS, NMEAS, NCON, OPTIONS )
%
%   Inputs:
%     SYS - an LTI system
%     NMEAS - dimension of the vector 'y' 
%     NCON  - dimension of the vector 'u' 
%     OPTIONS - structure with possible fields
%            .tol - tolerance for relative error (default 1.e-4)
%            .verbose - verbosity level (default 0)
%            .method  - solution method of
%                'gl' - formulas by Glover & Doyle (default)
%                'sa' - formulas by Safonov, Limebeer & Chiang 
%   Outputs:
%       K     - Hinf-optimal controller as an LTI model
%	    GAMMA - optimal value of gamma = ||Tzw||inf
%	    CL    - closed-loop system as an LTI model
%
%   See also H2REG, SDHINFREG, SDAHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check syntax
%------------------------------------------------------
    if ~exist('o2','var'), o2 = 1; end;
    if ~exist('i2','var'), i2 = 1; end;
    if exist('options','var'),
       if isfield(options,'tol'),     gammaTol = options.tol; end;
       if isfield(options,'verbose'), verbose = options.verbose; end;
       if isfield(options,'method'),  method = options.method; end;
    end;
    if ~exist('gammaTol', 'var'), gammaTol = 1.e-4; end;
    if ~exist('verbose', 'var'),  verbose = 0;      end;
    if ~exist('method', 'var'),   method = 'gl';    end;
    if ~isequal(method,'sa')  &&  ~isequal(method,'gl')
       error('Unknown solution method ''%s''', method); 
    end;    
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
    sys = minreal(ss(sys), sqrt(eps), 0);
    [valid,sys] = regular ( sys, o2, i2 );
   	%[a,b,c,d] = ssdata(sys);
   	%[n,n] = size(a);
  	[nout,nin] = size(sys);
	i1 = nin - i2;
	o1 = nout - o2;
%------------------------------------------------------
%   Define starting point	
%------------------------------------------------------
    sysH2 = sys;
    sysH2.d(1:o1,1:i1) = 0;
    Kh2 = h2reg ( sysH2, o2, i2 );    
    gamma = norm(lft(sys,Kh2), Inf);
    f = hinfstep ( sys, o1, o2, i2, gamma, method );
%------------------------------------------------------
%   Find bounds
%------------------------------------------------------
    if (f < 0)  &&  (f > -1) % gamma is admissible
       gamma1 = gamma; f1 = f;
      %------------------------------------------------
      %   Find lower bounds
      %------------------------------------------------
       while 1
         gamma = gamma / 2;
         f = hinfstep ( sys, o1, o2, i2, gamma, method );
         if (f < -1)  ||  (f > 0), break; end;
       end;
    else
       gamma1 = gamma;
      %------------------------------------------------
      %   Find upper bounds
      %------------------------------------------------
       while 1
         gamma1 = gamma1 * 2;
         f1 = hinfstep ( sys, o1, o2, i2, gamma1, method );
         if (f1 < 0)  &&  (f1 > -1), break; end;
       end;
    end;
%------------------------------------------------------
%   Bisection until positive definiteness of X and Y
%------------------------------------------------------
    while f < -1
      c = (gamma + gamma1) / 2;
      fc = hinfstep ( sys, o1, o2, i2, c, method );
      if (fc < 0)  &&  (fc > -1), 
        gamma1 = c; f1 = fc; 
      else
        gamma = c; f = fc;
        if fc > 0, break; end;
      end;
    end;
%--------------------------------------------------------------
%   Using gamma-iteration
%--------------------------------------------------------------
    ws = warning; warning off;
    c = gamma; count = 0;
    while abs(gamma - gamma1) > gammaTol
       c1 = c; 
       c = gamma + (gamma1 - gamma)*f/(f - f1);
       if abs(c-c1) < 0.02*(gamma1+gamma)
         c = (gamma + gamma1) / 2;
       end;
       if verbose,  fprintf('gamma = %g\n', gamma ); end;            
       fc = hinfstep ( sys, o1, o2, i2, c, method );
       [gamma gamma1 c fc], pause
       if (fc < 0)  &&  (fc > -1) 
         gamma1 = c; f1 = fc;
       else
         gamma = c; f = fc;
       end;       
       count = count + 1;
    end;        
    warning(ws);
    count
%--------------------------------------------------------------
%   Find optimal controller
%--------------------------------------------------------------    
    gamma = (gamma + gamma1)/2;
    while 1 
      [f,Kopt] = hinfstep ( sys, o1, o2, i2, gamma, method );
      if (f < 0)  &&  (f > -1), break; end;
      gamma = gamma + gammaTol/2;
    end;
%--------------------------------------------------------------
%   Constructing closed-loop system model
%--------------------------------------------------------------
    ClSys = minreal ( lft(sys, Kopt), sqrt(eps), 0 );

%------- End of HINFREG.M --------- KYuP ----------

%##############################################################
function [G,Kgamma] = hinfstep ( sys, o1, o2, i2, gamma, method )
%--------------------------------------------------------------
% HINFSTEP - find G = lam(P*S) - 1 for LTI system
%--------------------------------------------------------------    
%    Scale LTI model
%--------------------------------------------------------------
    sys.c(1:o1,:) = sys.c(1:o1,:) / gamma;
    sys.d(1:o1,:) = sys.d(1:o1,:) / gamma;
%------------------------------------------------------------
%    Try to construct gamma-suboptimal controller
%------------------------------------------------------------     
    ws = warning; %warning off    
    if isequal(method,'sa')
         [Kgamma,gamma0] = hinfone  ( sys, o2, i2 ); % --- Safonov, Limebeer & Chiang
    else [Kgamma,gamma0] = hinfone1 ( sys, o2, i2 ); % --- Glover & Doyle
    end;
    warning(ws);	
    G = gamma0 - 1;
       
%------- End of LTIXI.M --------- KYuP ----------
