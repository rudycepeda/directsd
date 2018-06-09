function [K,gamma] = hinfone ( sys, o2, i2 )
%HINFONE Find a subHinf-controller for an LTI-system with GAMMA=1.
%	     Formulas by Safonov, Limebeer and Chiang (Int. J. Control, 1989)	
%
%     [K,GAMMA] = HINFONE ( SYS, NMEAS, NCON )
%
%   Inputs:
%     SYS - an LTI system
%     NMEAS - dimension of the vector 'y' 
%     NCON  - dimension of the vector 'u' 
%
%   Outputs:
%     K     - Hinf-optimal controller as an LTI system such that
%             ||Tzw|| = 1
%     GAMMA = -1, if at least one Riccati solution is not positive definite
%           = max(lam(P*S)) otherwise
%

% References:
%  [1] M.G. Safonov, D.J.N. Limebeer, and R.Y. Chiang,  
%      Simplifying Hinf-theory via loop shifting, matrix 
%      pencil and descriptor concept, Int. J. Control, 1990, 
%      vol. 50, pp, 2467-2488.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 14-Jan-2006 $
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
    sys = minreal(ss(sys), sqrt(eps), 0);
    [valid,sys] = regular ( sys, o2, i2 );
   	[a,b,c,d] = ssdata(sys);
   	[nout,nin] = size(sys);
	i1 = nin - i2;
	o1 = nout - o2;
%------------------------------------------------------
%	Extract state-space matrices
%------------------------------------------------------
    [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
%------------------------------------------------------
%   Initialization
%------------------------------------------------------
    K = [];
%---------------------------------------------------------
%   If discrete system, use bilinear transformation
%---------------------------------------------------------
    if isdt(sys)
       sysb = bilinss ( sys );
       [K,gamma] = hinfone ( sysb, o2, i2 );
       if ~isempty(K)
         K    = bilinss ( K );        
         K.Ts = sys.Ts;
       end;
       return;     
    end;        
%---------------------------------------------------------
%   Scaling with respect D12 and D21
%---------------------------------------------------------
    d220 = d22;   
    d22 = d22 * 0;
    [b1,b2,c1,c2,d11,d12,d21,d22,su1,su2,sy1,sy2] = scaless (b1,b2,c1,c2,d11,d12,d21,d22);    
%---------------------------------------------------------
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
%   Find fopt such that det(I + fopt*sy2*d22*su2) = 0
%----------------------------------------------------------
    if i1 == o2 || o1 == i2
         fopt = - q22;
    else fopt  = -(q22 + q21 * inv(eye(i1-o2)-q11'*q11) * q11' * q12); 
    end;

    imax = 23;
    foptscale = 2^(-imax);
    fopt1 = fopt;
    for i = 1:imax+1
        f = eye(size(fopt,1)) + fopt1 * sy2 * d220 * su2;
        sigma = svd(f);
        if ~any(sigma < 1.e-8), break; end;
        sfopt = svd(fopt);
        fopt1 = (1 + foptscale*(1-gamma) / sfopt(1,1)) * fopt;
        foptscale = 2*foptscale;
    end
    fopt = fopt1;
    if ~(gamma < 1 && foptscale < 1.99)
%       error('Failed to find optimality point'); 
       gamma = -1;
       return;
    end
%----------------------------------------------------------
%   Zero out D11 term via loop shifting
%   (Safonov, Limebeer and Chiang, 1989, Int. J. Control)
%----------------------------------------------------------
    a   = a   + b2 * fopt * c2;
    b1  = b1  + b2 * fopt * d21;
    c1  = c1  + d12 * fopt * c2;
    d11 = d11 + d12 * fopt * d21;

    X = d11;
    IXX1 = eye(o1)- X * X';
    IXX2 = eye(i1)- X' * X;
    IIXX1 = inv(IXX1);

    a   = a + b1 * X' * IIXX1 * c1;
    b2  = b2 + b1 * X' * IIXX1 * d12;
    b1  = b1 * IXX2^(-0.5);
    c2  = c2 + d21 * X' * IIXX1 * c1;
    c1  = IXX1^(-0.5) * c1;
    d22 = d21 * X' * IIXX1 * d12;

    d11 = d11 * 0;
    d12 = IXX1^(-0.5) * d12;
    d21 = d21 * IXX2^(-0.5);  
%-------------------------------------------------------- 
%   Zero out the new D22
%----------------------------------------------------------
    d22a = d22;
    d22 = d22 * 0;
%---------------------------------------------------------
%   Scaling again with respect to new D12 and D21
%----------------------------------------------------------
    [b1,b2,c1,c2,d11,d12,d21,d22,tu1,tu2,ty1,ty2] = scaless (b1,b2,c1,c2,d11,d12,d21,d22);    
%----------------------------------------------------------
%	Solution to the state feedback problem (P-Riccati)
%----------------------------------------------------------
    c1til = (eye(o1) - d12*d12')*c1;
    b1til = b1 * (eye(i1) - d21'*d21);

    a1 = a - b2*d12'*c1;
    eigA1 = eig(a1);
%----------------------------------------------------------
%   Special case: P = 0, if the Hamiltonian A matrix is stable 
%   and D12 is square
%----------------------------------------------------------
    if (o1 == i2)  &&  all(real(eigA1)< 0)
       P2 = zeros(n,n);
       P1 = eye(n);
       P  = zeros(n,n);
       eigP_inf = eigA1;
%----------------------------------------------------------
%   General case
%----------------------------------------------------------
    else     
        q1 = c1til'*c1til;
        r1 = -(b1*b1' - b2*b2');
        [P,Poles,err,P1,P2] = care2 ( a1, q1, r1 );
        [vinf,dinf] = eig ( a1*P1 - b2*b2'*P2, P1 );
        eigP_inf = diag(dinf);
        ind = find ( isinf(eigP_inf) | isnan(eigP_inf) );
        eigP_inf(ind) = 1/eps;        
        if max(real(eigP_inf)) > 0
%           error('Matrix P is not positive definite'); 
           gamma = -1;
           return;
        end;                   
    end
%------------------------------------------------------
%	Solution to the output injection problem (S-Riccati)
%----------------------------------------------------------
    a2 = a - b1*d21'*c2;
    eigA2 = eig(a2);   
%------------------------------------------------------
%   Special case: S = 0, if the Hamiltonian A matrix is stable 
%   and D21 is square
%------------------------------------------------------
    if (o2 == i1)  &&  all(real(eigA2) < 0)
       S2 = zeros(n,n);
       S1 = eye(n); 
       S  = zeros(n,n);
       eigS_inf = eigA2;
%------------------------------------------------------
%   General case
%------------------------------------------------------
    else 
       q2  = b1til*b1til';
       r2  = -(c1'*c1 - c2'*c2);
       [S,Poles,err,S1,S2] = care2 ( a2', q2, r2 );
       S1 = S1';
       S2 = S2';       
       [vinf,dinf] = eig ( S1*a2 - S2*c2'*c2, S1 );
       eigS_inf = diag(dinf);
       ind = find ( isinf(eigS_inf) | isnan(eigS_inf) );
       eigS_inf(ind) = 1/eps;
       if max(real(eigS_inf)) > 0
%           error('Matrix S is not positive definite'); 
           gamma = -1;
           return;
       end;                    
    end     
%------------------------------------------------------
%	Check solution
%------------------------------------------------------
    [v,dd] = eig (P2'*S2', P1'*S1');
    eigPS = diag(dd);
    ind = find ( isinf(eigPS) | isnan(eigPS) );
    eigPS(ind) = 1/eps;
    gamma = max(real(eigPS));
    if gamma > 1
%       error('No stabilizing H-infinity controller exists');  
       return;
    end;
%-------------------------------------------------------- 
%   Controller parameterization K(s) in descriptor form
%------------------------------------------------------
    E = S1*P1-S2*P2;

    ak  = a-b2*d12'*c1-b1*d21'*c2;
    ak  = S1*ak*P1+S2*ak'*P2 + S1*(b1til*b1til'-b2*b2')*P2 + ...
          S2*(c1til'*c1til-c2'*c2)*P1;
       
    bk1 = S2*c2' + S1*b1*d21';
    bk2 = S1*b2 + S2*c1'*d12;
    ck1 = -(b2'*P2+d12'*c1*P1);
    ck2 = -(c2*P1 + d21*b1'*P2);

    ki1 = size(bk1, 2);
    ki2 = size(bk2, 2);
    ko1 = size(ck1, 1);
    ko2 = size(ck2, 1);

    dk11 = zeros(ko1,ki1);
    dk12 = eye(ko1);
    dk21 = eye(ko2);
    dk22 = zeros(ko2,ki2);    
%------------------------------------------------- 
%   Reverse the controller scaling (stage II)
%------------------------------------------------------
    bk1  = bk1*ty2;
    ck1  = tu2*ck1;
    dk11 = tu2*dk11*ty2;
    dk12 = tu2*dk12;
    dk21 = dk21*ty2;
%-------------------------------------------------- 
%   Shifting D22_A
%------------------------------------------------------
    temp = inv( eye(size(dk11,1)) + dk11*d22a );
    ak   = ak - bk1 * d22a * temp * ck1;
    bk2  = bk2 - bk1 * d22a * temp * dk12;
    ck2  = ck2 - dk21 * d22a * temp * ck1;
    ck1  = temp * ck1;
    dk22 = dk22 - dk21 * d22a * temp * dk12;
    dk12 = temp * dk12;
    dk11 = temp * dk11;

    temp = eye (size(d22a,1)) - d22a * dk11;
    bk1  = bk1 * temp;
    dk21 = dk21 * temp;
%--------------------------------------------------- 
%   Reverse the "fopt" term:
%------------------------------------------------------
    dk11 = dk11 + fopt;
%---------------------------------------------------
%   Reverse the controller scaling (stage I)
%------------------------------------------------------
    bk1  = bk1 * sy2;
    ck1  = su2 * ck1;
    dk11 = su2 * dk11 * sy2;
    dk12 = su2 * dk12;
    dk21 = dk21 * sy2;
%--------------------------------------------------- 
%   Shifting initial D22 term
%------------------------------------------------------
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
%--------------------------------------------------------  
    if cond(E) > 1/eps
%      error('Singular matrix E in descriptor form');        
       gamma = 1; % --- singular solution
    end;        
    [ak,bk,ck,dk,OK] = descr2ss ( ak, [bk1, bk2], [ck1; ck2], [dk11 dk12; dk21 dk22], E );
    if ~OK
       warning('Controller is improper: erroneous results possible'); 
    end;        
                
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
%------------------------------------------------------
    ac = ak; bc = bk1; cc = ck1; dc = dk11;
    K = ss ( ac, bc, cc, dc );		
    K = minreal ( K, sqrt(eps), 0 );
    
%------- End of HINFONE.M --------- KYuP ----------


