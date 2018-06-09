function [K,gamma] = hinfone1 ( sys, o2, i2, gamma )
%HINFONE1 Fund subHinf-controller for an LTI-system.
%         Formulas by K. Glover  &  J. Doyle  (1988).
%
%     [K,GAMMA] = HINFONE1 ( SYS, NMEAS, NCON, GAMMA )
%
%   Inputs:
%     SYS - an LTI system
%     NMEAS - dimension of the vector 'y' 
%     NCON  - dimension of the vector 'u' 
%     GAMMA - gamma (default 1)
%
%   Outputs:
%     K    - Hinf-optimal controller as an LTI system such that
%            ||Tzw|| <= gamma, or empty matrix if none exists
%     GAMMA - success indicator (0: OK; -1: no solution)
%

% References:
%  [1] K. Glover, and J.C. Doyle, "State-space formulae for all 
%      stabilizing controllers that satisfy an Hinf-norm bound 
%      and relations to risk sensitivity," Systems and Control 
%      Letters, 1988. vol. 11, pp. 167-172.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%   Use bilinear transformation for a discrete-time system 
%------------------------------------------------------
    if isdt(sys)
       sysb = bilinss ( sys );
       [K,gamma] = hinfone1 ( sysb, o2, i2 );
       if ~isempty(K)
         K = bilinss ( K );        
         K.Ts = sys.Ts;
       end;
       return;     
    end;        
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
	sys = minreal(ss(sys), sqrt(eps), 0);
    [valid,sys] = regular ( sys, o2, i2 );
    [a,b,c,d] = ssdata(sys);
    [n,n] = size(a);
    [nout,nin] = size(sys);
	i1 = nin - i2;
	o1 = nout - o2;
%------------------------------------------------------
%	Extract state-space matrices
%------------------------------------------------------
    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
%------------------------------------------------------
%	Scaling
%------------------------------------------------------
    if ~exist('gamma','var'), gamma = 1; end;
    c1 = c1 / gamma;  
    d11 = d11 / gamma;  
    d21 = d21 / gamma;  
%------------------------------------------------------
%	Check by means of Robust Toolbox
%------------------------------------------------------
	if 0  &&  exist('mksys.m','file')	   
	   if isct(sys)
           save;
           A = a; B1 = b1; B2 = b2; C1 = c1; C2 = c2; 
           D11 = d11; D12 = d12; D21 = d21; D22 = d22;       
           hinfkgjd; 
         try  
           ac = acp; bc = bcp; cc = ccp; dc = dcp;
           disp(' ');
	       disp('Optimal controller as given by Robust Toolbox:');
           disp('----------------------------------------------');            
	       C0 = zpk(ss(ac,bc,cc,dc,sys.Ts))
           dcgain(C0)
         catch end;  
           clear all;
           load;
       end;       
   end;
   
%------------------------------------------------------
%	Check by means of Mu-Tools Toolbox
%------------------------------------------------------
	if 0  &&  isct(sys)  &&  exist('hinfsyn.m','file')   
       plant = pck (sys.a, sys.b, sys.c, [d11 d12; d21 d22] );
       gmin = 0;
       gmax = 100;
       [k,g,gfin,ax,ay,hamx,hamy] = hinfsyn(plant, 1, 1, gmin, gmax, 1.e-6, 2, 1.e-10, 1.e-6, 0);
      try 
       disp(' ')
       disp('Optimal controller as given by Mu-Toolbox:');
       disp('-------------------------------------------');
	   [ac,bc,cc,dc] = unpck (k);
	   zpk(ss(ac,bc,cc,dc,sys.Ts))
      catch end; 
    end;   
%---------------------------------------------------------
%   Scaling with respect D12 and D21
%---------------------------------------------------------
    d220 = d22; 
    d22 = d22 * 0;
    [b1,b2,c1,c2,d11,d12,d21,d22,su1,su2,sy1,sy2] = scaless (b1,b2,c1,c2,d11,d12,d21,d22);    
%------------------------------------------------------
%	Initialization
%------------------------------------------------------
    K = [];
    gamma = [];
% ---------------------------------------------------------
%   Find lower bound for Hinf norm
%---------------------------------------------------------
    q11 = d11(1:o1-i2,1:i1-o2);
    q21 = d11(o1-i2+1:o1,1:i1-o2);
    q12 = d11(1:o1-i2,i1-o2+1:i1);
    q22 = d11(o1-i2+1:o1,i1-o2+1:i1);    
    gammar = max([svd([q11 q12]);0]);
    gammac = max([svd([q11;q21]);0]);
    gamma = max(gammar, gammac);    
    if gamma > 1, return; end;          % no solution
%----------------------------------------------------------
%   Start computations: use general formulas with gamma=1    
%---------------------------------------------------------
    gam = 1;
    gam2 = gam*gam;
%---------------------------------------------------------- 
%   Fixing D11 matrix
%---------------------------------------------------------
    if i1 == o2  ||  o1 == i2
       dk11 = - q22;
       dk12 = eye(i2);
       dk21 = eye(o2);
    else
       gamd1  = gam2*eye(o1-i2) - q11*q11';
       gamd2  = gam2*eye(i1-o2) - q11'*q11;
       dk11   = - q21*q11'*inv(gamd1)*q12 - q22;
       dk1212 = eye(i2) - q21*inv(gamd2)*q21';
       dk12   = chol(dk1212); dk12 = dk12';
       dk2121 = eye(o2) - q12'*inv(gamd1)*q12;
       dk21   = chol(dk2121);
    end
%----------------------------------------------------------
%   Solving X Riccati equation
%--------------------------------------------------------- 
    Zr = zeros(i1+i2);
    Zr(1:i1,1:i1) = gam2*eye(i1);
    d1d = [d11 d12];
    R = d1d'*d1d - Zr;
    Ri = inv(R);
    B = [b1 b2];
    ax = a - B*Ri*d1d'*c1;
    rx = B*Ri*B';
    qx = c1'*c1 - c1'*d1d*Ri*d1d'*c1;
    
    [X,Poles,err,X1,X2] = care2 ( ax, qx, rx );
%----------------------------------------------------------
%   Solving Y Riccati equation
%---------------------------------------------------------  
    Zrw = zeros(o1+o2, o1+o2);
    Zrw(1:o1,1:o1) = gam2*eye(o1);
    dd1 = [d11; d21];
    Rw = dd1*dd1' - Zrw;
    Rwi = inv(Rw);
    C = [c1; c2];
    ay = a' - C'*Rwi*dd1*b1';
    ry = C'*Rwi*C;
    qy = b1*b1' - b1*dd1'*Rwi*dd1*b1';

    [Y,Poles,err,Y1,Y2] = care2 ( ay, qy, ry );
%----------------------------------------------------------
%   Check solvability:
%      1) matrices X and Y must be positive definite
%      2) ro(X*Y) < 1 
%---------------------------------------------------------
    lamTol = 1e-8;
    lamX = min(real(eig(X)));
    lamY = min(real(eig(Y)));
    gamma = max(real(eig(X*Y)));
    if gamma < 0, gamma = 1 + abs(gamma); end;
    if (lamX<0 && abs(lamX)>lamTol)  || ... % X and/or Y is not positive definite 
       (lamY<0 && abs(lamY)>lamTol)
      gamma = -1; 
      return;
    end;
    if gamma > 1
%       error('No stabilizing controller exists'); 
       return;  
    end;
%---------------------------------------------------------- 
%   Auxiliary matrices
%---------------------------------------------------------
    F   = - Ri*(d1d'*c1+B'*X);
    F11 = F ( 1:(i1-o2), : );
    F12 = F ( (i1-o2+1):i1, : );
    F2  = F ( (i1+1):(i1+i2),:);

    H   = - (b1*dd1'+Y*C')*Rwi;
    H11 = H ( :, 1:(o1-i2) );
    H12 = H ( :, (o1-i2+1):o1 );
    H2  = H ( :, (o1+1):(o1+o2) );
%----------------------------------------------------------- 
%   Hinf-optimal controller in descriptor form
%---------------------------------------------------------
    Z   = eye(n) - inv(gam2)*Y*X;
    bk2 = (b2 + H12)*dk12;
    ck2 = - dk21*(c2 + F12);
    bk1 = - H2 + bk2*inv(dk12)*dk11;
    ck1 = F2 + dk11*inv(dk21)*ck2;
    ak  = a*Z + H*C*Z + bk2*inv(dk12)*ck1;
%---------------------------------------------------
%   Reverse the controller scaling (stage I)
%---------------------------------------------------------
    bk1  = bk1 * sy2;
    ck1  = su2 * ck1;
    dk11 = su2 * dk11 * sy2;
    dk12 = su2 * dk12;
    dk21 = dk21 * sy2;
%--------------------------------------------------- 
%   Shifting initial D22 term
%---------------------------------------------------------
    dk22 = zeros( size(dk21,1), size(dk12,2) );
    temp = inv( eye(size(dk11,1)) + dk11 * d220 );
    ak   = ak - bk1 * d220 * temp * ck1;
    bk2  = bk2 - bk1 * d220 * temp * dk12;
    ck2  = ck2 - dk21 * d220 * temp * ck1;
    ck1  = temp * ck1;
    dk22 = dk22 - dk21 * d220 * temp * dk12;
    dk12 = temp * dk12;
    dk11 = temp * dk11;

    temp = eye(size(d220,1)) - d220*dk11;
    bk1  = bk1 * temp;
    dk21 = dk21 * temp;
%--------------------------------------------------------  
%   Putting descriptor controller K(s) in state space form:
%---------------------------------------------------------
    if cond(Z) > 1/eps
%      error('Singular matrix E in descriptor form');                
       return;
    end;    
    [ak,bk,ck,dk,OK] = descr2ss ( ak, [bk1, bk2], [ck1; ck2], [dk11 dk12; dk21 dk22], Z );
    if ~OK
       warning('Controller is improper: erroneous results possible'); 
    end;        
                
    ki1 = size(bk1, 2);
    ki2 = size(bk2, 2);
    ko1 = size(ck1, 1);
    ko2 = size(ck2, 1);

    %ki = ki1 + ki2;                
    %ko = ko1 + ko2;                
    
    bk1  = bk ( :, 1:ki1 );
    %bk2  = bk ( :, (ki1+1):ki );
    ck1  = ck ( 1:ko1, : );
    %ck2  = ck ( (ko1+1):ko, : );
    dk11 = dk ( 1:ko1, 1:ki1 );
    %dk12 = dk ( 1:ko1, (ki1+1):ki );
    %dk21 = dk ( (ko1+1):ko, 1:ki1 );
    %dk22 = dk ( (ko1+1):ko, (ki1+1):ki );
%------------------------------------------------------
%	Final Hinf-optimal controller
%---------------------------------------------------------
    ac = ak; bc = bk1; cc = ck1; dc = dk11;
    K = ss ( ac, bc, cc, dc );		
	K = minreal(K, sqrt(eps), 0);
    
%------- End of HINFONE1.M --------- KYuP ----------
