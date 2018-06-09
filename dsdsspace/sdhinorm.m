function [gamma,gamD11,count,gamDT] = sdhinorm ( csys, ksys, tol, type )
%SDHINORM Hinf-norm (L2-induced norm) for a sampled-data system.
%
%     [GAMMA,GAMD11,COUNT,GAMDT] = SDHINORM ( CSYS, KSYS, TOL )
%     [GAMMA,GAMD11,COUNT,GAMDT] = SDHINORM ( CSYS, KSYS, METHOD )
%     [GAMMA,GAMD11,COUNT,GAMDT] = SDHINORM ( CSYS, KSYS, TOL, METHOD )
%
%   Inputs:
%     CSYS - continuous-time standard system 
%	  KSYS - discrete-time controller (positive feedback)
%     TOL  - tolerance (default 1e-4)
%	  METHOD - type of formulas for Hinf-equivalent model: 
%         'ch' - Chen and Francis
%         'ba' - Bamieh and Pearson
%         'ca' - Cantoni and Glover
%         'ha' - Hayakawa, Hara and Yamamoto
%         'mi' - Mirkin and Palmor (default)
%   Outputs:
%	  GAMMA  - Hinf-norm of the closed-loop system
%	  GAMD11 - Hinf-norm of the operator D11
%     COUNT  - number of iteration steps 
%	  GAMDT  - Hinf-norm of the discretized system
%
%   See also SDHINFREG, SDAHINORM, SDH2NORM.

% References:
% [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%     Springer-Verlag, Berlin etc., 1995.
% [2] Bamieh B.A., Pearson J.B., A general framework for linear periodic systems
%     with applications to Hinf sampled-data control // IEEE
%     Trans. Automat. Contr., vol. AC-37, no. 4, pp. 418-435, 1992.
% [3] M. W. Cantoni, K. Glover, Hinf sampled-data synthesis and 
%     related numerical issues, Automatica, vol. 33, no. 12, 
%     pp. 2233-2241, 1997.  
% [4] Hayakawa Y., Hara S., and Yamamoto Y., Hinf type problem for
%     sampled-data control systems - a solution via minimum energy 
%     characterization // IEEE Trans. Automat. Contr., vol. AC-39, no. 11, 
%     pp. 2278-2284, 1994.
% [5] Mirkin L., Tadmor G., Yet another Hinf-discretization //
%     IEEE Trans. Automat. Contr., vol. AC-48, no. 5, pp. 891--894, 2003.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check input data
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1e-4; 
        elseif ischar(tol)
          type = tol; tol = 1e-6;  
        end;
        if ~exist('type','var'), type = 'mi'; end; 
        tolmin = sqrt(eps);
        csys = minreal(ss(csys),tolmin,0);
        [oc,ic] = size ( csys );
        if isempty(ksys),
           gamma = norm(csys, 'inf');
           gamD11 = NaN;
           count = 0;
           return;
        else
    	   if ~isdt(ksys)
    	     error('Controller must be a discrete time model');
    	   end;
           ksys = minreal(ss(ksys),tolmin,0);
           T = ksys.Ts;
           [ok,ik] = size ( ksys );
    	   o1 = oc - ik;
    	   i1 = ic - ok;
        end;        
    	if ~isct(csys)
    	   error('Generalized plant must be a continuous time model');
    	end;
    	if i1 < 1  ||  o1 < 1
    	   error('Incorrect plant-controller pair');
    	end;
    	b1 = csys.b ( :, 1:i1 );
    	c1 = csys.c ( 1:o1, : );
    	d1 = csys.d ( 1:o1, 1:i1 );
%------------------------------------------------------
%	Find norm of the operator D11	
%------------------------------------------------------
%        if nargout > 1,
%          gamD11 = d11norm ( csys.a, b1, c1, d1, T );
%        end;
%------------------------------------------------------
%	Initial guess - Hinf norm of the discretized system
%------------------------------------------------------
        %dsys = c2d ( csys, T );
        Nfast = 5;
        [i2,o2] = size( ksys );
        dsys = sdfast ( csys, T, Nfast, o2, i2, 'l2' );
        clFast = lft ( dsys, ksys );
        gamma = norm ( clFast, 'inf');        
        f = normMinus1 ( csys, ksys, gamma, type );
%------------------------------------------------------
%       Search for an uper bound
%------------------------------------------------------
        gamma1 = gamma;
        while 1
          gamma1 = 2*gamma1;
          f1 = normMinus1 ( csys, ksys, gamma1, type );
          if f1 < 0, break; end;
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
           count = count + 1;  
           fc = normMinus1 ( csys, ksys, c, type ); 
           if (fc < 0)  &&  (fc > -1) 
             gamma1 = c; f1 = fc;
           else
             gamma = c; f = fc;
           end;
        end;

%------- End of SDHINORM.M --------- KYuP ----------

function OK = checkone ( csys, ksys, gamma, type )
%CHECKONE Check if Hinf-norm of SD system is less than GAMMA.
%
%     OK = CHECKONE ( CSYS, KSYS, GAMMA, TYPE )
%
%   Inputs:
%     CSYS  - continuous-time plant
%     KSYS  - discrete-time controller
%     GAMMA - attenuation level
%	  TYPE  - type of formulas:
%         'ch' - Chen and Francis
%         'ba' - Bamieh and Pearson
%         'ca' - Cantoni and Glover
%         'ha' - Hayakawa, Hara and Yamamoto
%         'mi' - Mirkin and Palmor (default)
%
%   Outputs:
%      OK - is Hinf-norm of SD-system less than GAMMA (1/0)
%

% References:
%  [1] T. Chen, B.A. Francis, Optimal Sampled-Data Control Systems, 
%      Berlin etc.: Springer-Verlag, 1995.
%  [2] Mirkin L., Tadmor G., Yet another Hinf-discretization //
%      IEEE Trans. Automat. Contr., vol. AC-48, no. 5, pp. 891--894, 2003.
%  [3] Modification of the basic idea by K. Polyakov, (C) 2005
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        OK = 0;
%------------------------------------------------------
%   Construct Hinf-equivalent discrete-time model
%------------------------------------------------------
        [i2,o2] = size ( ksys );
        dsys = sdhinfmod ( csys, ksys.Ts, o2, i2, gamma, type );
        gamma = 1;
%------------------------------------------------------
%   Construct closed-loop system matrices
%------------------------------------------------------
        clSys = lft ( dsys, ksys );
        [Ad,Bd,Cd,Dd] = ssdata ( clSys );        
%------------------------------------------------------
%   Using bilinear transform
%------------------------------------------------------
        In = eye(size(Ad));
        On = In*0;
        Db = Cd/(Ad+In)*Bd;
        dGamma = norm(Db);
        if dGamma > gamma, return; end; % OK = 0
%------------------------------------------------------
%   Check eigs of matrix pencil
%------------------------------------------------------        
        [p,m] = size(Dd);
        R = gamma^2*eye(p) - Dd*Dd';
        S = gamma^2*eye(m) - Dd'*Dd;
        Adm = Ad + Bd*Dd'/R*Cd;
        CC = Cd'/R*Cd;
        BB = Bd/S*Bd';
        M = [Adm On; -CC*gamma In];
        N = [In -BB*gamma; On Adm'];
        eigMN = eig(M, N);
        ucEigs = find(abs(abs(eigMN)-1) < 1e-6);
        OK = isempty(ucEigs);

%------- End of CHECKONE.M --------- KYuP ----------

function G = normMinus1 ( csys, ksys, gamma, type )
%NORMMINUS1 Check if Hinf-norm of SD system is less than gamma.
%
%     G = NORMMINUS1 ( CSYS, KSYS, GAMMA, TYPE )
%
%   Inputs:
%     CSYS  - continuous-time plant
%     KSYS  - discrete-time controller
%     GAMMA - attenuation level
%	  TYPE  - type of formulas:
%         'ch' - Chen and Francis
%         'ba' - Bamieh and Pearson
%         'ca' - Cantoni and Glover
%         'ha' - Hayakawa, Hara and Yamamoto
%         'mi' - Mirkin and Palmor (default)
%
%   Outputs:
%      G - Hinf-norm of equivalent D-system minus 1
%

% References:
%  [1] T. Chen, B.A. Francis, Optimal Sampled-Data Control Systems, 
%      Berlin etc.: Springer-Verlag, 1995.
%  [2] Mirkin L., Tadmor G., Yet another Hinf-discretization //
%      IEEE Trans. Automat. Contr., vol. AC-48, no. 5, pp. 891--894, 2003.
%  [3] Modification of the basic idea by K. Polyakov, (C) 2005
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%   Construct Hinf-equivalent discrete-time model
%------------------------------------------------------
        [i2,o2] = size ( ksys );
        dsys = sdhinfmod ( csys, ksys.Ts, o2, i2, gamma, type );
%------------------------------------------------------
%   Construct closed-loop system matrices
%------------------------------------------------------
        clSys = lft ( dsys, ksys );
        [Ad,Bd,Cd,Dd] = ssdata ( clSys );        
        G = norm ( clSys, 'inf' ) - 1;

%------- End of NORMMINUS1.M --------- KYuP ----------

