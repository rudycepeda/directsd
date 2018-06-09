function [Kopt,gamma,Poles] = sdhinfreg ( csys, T, o2, i2, options )
%SDHINFREG Hinf-optimal controller for sampled-data system.
%
%     [K,GAMMA,POLES] = SDHINFREG ( SYS, T, NMEAS, NCON, OPTIONS )
%
%   Inputs:
%     SYS - an LTI system
%     T   - sampling period
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1)
%     OPTIONS - structure with possible fields
%          .tol - tolerance for relative error (default 1.e-4)
%          .verbose - verbosity level (default 0)
%          .type - type of discrete-time equivalent model:
%              'ch' - Chen and Francis
%              'ba' - Bamieh and Pearson
%              'ca' - Cantoni and Glover
%              'ha' - Hayakawa, Hara and Yamamoto
%              'mi' - Mirkin and Palmor (default)
%          .method - iteration method:
%              'sa' - formulas by Safonov, Limebeer & Chiang 
%              'gl' - formulas by Glover & Doyle (default)
%   Outputs:
%     K     - Hinf-optimal controller as an LTI model
%	  GAMMA - optimal value of gamma = ||Tzw||inf
%	  POLES - poles of the closed-loop system in 'z'-plane
%
%   See also SDH2REG, SDH2, SDAHINF, SDHINORM.

% References:
% [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%     Springer-Verlag, Berlin etc., 1995.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check syntax
%------------------------------------------------------
    if exist('options','var'),
       if isfield(options,'tol'),     tol = options.tol; end;
       if isfield(options,'verbose'), verbose = options.verbose; end;
       if isfield(options,'type'),    type = options.type; end;
       if isfield(options,'method'),  method = options.method; end;
    end;
    if ~exist('tol', 'var'),      tol = 1.e-4; end;
    if ~exist('verbose', 'var'),  verbose = 0;      end;
    if ~exist('type', 'var'),     type = 'mi';      end;
    if ~exist('method', 'var'),   method = 'sa';    end;
    if ~isequal(method,'sa')  &&  ~isequal(method,'gl')
       error('Unknown solution method ''%s''', method); 
    end;    
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
    csys = minreal(ss(csys), sqrt(eps), 0);
  	[nout,nin] = size(csys);    
    if ~exist('i2','var'), i2 = 1; end;
    if ~exist('o2','var'), o2 = 1; end;
    [o1,i1] = checkio ( nout, nin, o2, i2 );
%------------------------------------------------------
%	Find norm of the operator D11	
%------------------------------------------------------
  	b1 = csys.b ( :, 1:i1 );
   	c1 = csys.c ( 1:o1, : );
   	d1 = csys.d ( 1:o1, 1:i1 );
    gamD11 = d11norm ( csys.a, b1, c1, d1, T );
%------------------------------------------------------
%   Initial guess: fast sampled system	
%------------------------------------------------------
    %N = 10;
    %dsys = sdfast ( csys, T, N, o2, i2 );
    %[Kd,gamDT] = hinfreg ( dsys );
    gamDT = gamD11;
%------------------------------------------------------
%	Find bounds using bisection
%------------------------------------------------------
    gamma = gamDT;
    f = sdhinfstep ( csys, T, o2, i2, gamma, type, method );
    gamma1 = gamma;
    while 1
      gamma1 = 2*gamma1;
      f1 = sdhinfstep ( csys, T, o2, i2, gamma1, type, method );
      if (f1 < 0)  &&  (f1 > -1), break; end;
    end;
%------------------------------------------------------
%   Bisection until positive definiteness of X and Y
%------------------------------------------------------
    while 1
      c = (gamma + gamma1) / 2;
      fc = sdhinfstep ( csys, T, o2, i2, c, type, method );
      if (fc < 0)  &&  (fc > -1), 
        gamma1 = c; f1 = fc; 
      else
        gamma = c; f = fc;
        if fc > 0, break; end;
      end;
    end;
%------------------------------------------------------
%       Main iteration loop 
%------------------------------------------------------
    count = 0;
    c = gamma;
    while abs(gamma - gamma1) > tol
       c1 = c; 
       c = gamma + (gamma1 - gamma)*f/(f - f1);
       if abs(c-c1) < 0.02*(gamma1+gamma)
         c = (gamma + gamma1) / 2;
       end;
       if verbose,  fprintf('gamma = %g\n', gamma ); end;            
       count = count + 1;  
       fc = sdhinfstep ( csys, T, o2, i2, c, type, method );
       if (fc < 0)  &&  (fc > -1) 
         gamma1 = c; f1 = fc;
       else
         gamma = c; f = fc;
       end;
    end;
%--------------------------------------------------------------
%   Find optimal controller
%--------------------------------------------------------------
    gamma = (gamma + gamma1)/2;
    while 1 
      [f,Kopt] = sdhinfstep ( csys, T, o2, i2, gamma, type, method ); 
      if (f < 0)  &&  (f > -1), break; end;
      gamma = gamma + tol/2;
    end;
%--------------------------------------------------------------
%   Constructing closed-loop system model
%--------------------------------------------------------------
    if i2 == 1  &&  o2 == 1,
      Kopt = ss(minreal(zpk(Kopt), 1e-3));
    end;
    Pd = c2d ( csys, T );
    clSys = minreal ( lft ( Pd, Kopt ), sqrt(eps), 0 );
    Poles = eig ( clSys.a );

%------- End of SDHINFREG.M --------- KYuP ----------

%##############################################################
function [G,Kgamma] = sdhinfstep ( csys, T, o2, i2, gamma, type, method )
%SDHINFSTEP Design Hinf-equivalent discrete model and
%           Hinf-optimal discrete controller
%           
%--------------------------------------------------------------
%    Construct equivalent discrete model for a current gamma
%--------------------------------------------------------------
    dsys = sdhinfmod ( csys, T, o2, i2, gamma, type );
%------------------------------------------------------------
%    Try to construct gamma-suboptimal controller
%------------------------------------------------------------     
    if isequal(method,'sa')
         [Kgamma,gamma0] = hinfone  ( dsys, o2, i2 ); % --- Safonov, Limebeer & Chiang
    else [Kgamma,gamma0] = hinfone1 ( dsys, o2, i2 ); % --- Glover & Doyle
    end;
    G = gamma0 - 1;
       
%------- End of SDHINFSTEP.M --------- KYuP ----------
