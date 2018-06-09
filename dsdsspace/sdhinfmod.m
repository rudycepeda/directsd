function [dsys,DP11,DP12,DP21,DP22] = sdhinfmod ( csys, T, o2, i2, gamma, type )
%SDHINFMOD Hinf-discrete equivalent of a sampled-data system.
%
%     [DSYS,DP11,DP12,DP21,DP22] = SDHINFMOD ( SYS, T, NMEAS, NCON, GAMMA, TYPE )
%
%   Inputs:
%     SYS - continuous-time standard system 
%               | z | = | P11 P12 | | w |
%               | y |   | P21 P22 | | u |
%  	  T     - sampling period	
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1)
%	  GAMMA - bound on induced Hinf-norm  (default 1)
%	  TYPE  - type of formulas: 
%        'ch' - Chen and Francis
%        'ba' - Bamieh and Pearson
%        'ca' - Cantoni and Glover
%        'ha' - Hayakawa, Hara and Yamamoto
%        'mi' - Mirkin and Palmor (default)
%
%   Outputs:
%	  DSYS  - discrete-time Hinf-equivalent of the system
%               | zd | = | DP11 DP12 | | dw |
%               | yd |   | DP21 DP22 | | du |
%

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
%	Check type of model
%------------------------------------------------------
    if ~exist('type','var'), type = 'mi'; end;
    if ~isequal(type,'ba')  &&  ~isequal(type,'ca') ...
       &&  ~isequal(type,'mi')  &&  ~isequal(type,'ch')  &&  ~isequal(type,'ha')         
       error('Unknown type of mode ''%s''', type); 
    end;
%------------------------------------------------------
%	Continuous-time state space
%------------------------------------------------------
	sys = minreal(ss(csys), sqrt(eps), 0);
	if ~isct(sys), error('Initial system must be a continuous time model'); end;
    [a,b,c,d] = ssdata(sys);
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
    
    if ~exist('gamma','var'), gamma = 1; end;
    gamma2 = gamma^2;
    xi = 1/gamma; xi2 = 1/gamma^2;
%------------------------------------------------------
%	Extract partial matrices and check assumptions
%------------------------------------------------------
	[b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( b, c, d, o2, i2 );
	if norm(d21) > eps, warning('D21 must be (assumed) zero'); end;
	if norm(d22) > eps, warning('D22 must be (assumed) zero'); end;
	%d21 = 0*d21;
	%d22 = 0*d22;       
%------------------------------------------------------
%	Remove d11 by loop shifting
%------------------------------------------------------
	if norm(d11) > eps, 
       % warning('D11 must be (assumed) zero'); 
       Rxi = eye(i1) - xi2*d11'*d11; sqrtRxi = sqrtm(Rxi);
       Sxi = eye(o1) - xi2*d11*d11'; sqrtSxi = sqrtm(Sxi);
       a = a + xi2*b1*d11'*(Sxi\c1);
       b1 = b1/sqrtRxi;
       b2 = b2 + xi2*b1*d11'*(Sxi\d12);
       c1 = sqrtSxi \ c1;
       d12 = sqrtSxi \ d12;
    end;
%------------------------------------------------------
%	Minimal realization of required blocks
%------------------------------------------------------  
    switch type
      case 'ch', 
        [Ad,Bd,Cd,Dd] = ChenFrancis ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 );   
      case 'ba', 
        [Ad,Bd,Cd,Dd] = BamiehPearson ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 );   
      case 'ca', 
        [Ad,Bd,Cd,Dd] = CantoniGlover ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 );   
      case 'mi', 
        [Ad,Bd,Cd,Dd] = MirkinPalmor ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 );   
      case 'ha', 
        [Ad,Bd,Cd,Dd] = HayakawaHaraYamamoto ( a, b1, b2, c1, c2, d12, T, xi, i2, o2 ); % xi instead of xi2           
    end;
    try
       dsys = minreal ( ss ( Ad, Bd, Cd, Dd, T ), sqrt(eps), 0 );
    catch keyboard; 
    end;        
    if nargout > 1,
      DP11 = minreal ( ss( Ad, B1d, C1d, zeros(o1d,i1d), T ), sqrt(eps), 0);
    end;
    if nargout > 2,
      DP12 = minreal ( ss( Ad, B2d, C1d, D12d, T ), sqrt(eps), 0);
    end;
    if nargout > 3,
      DP21 = minreal ( ss( Ad, B1d, c2, zeros(o2,i1d), T ), sqrt(eps), 0);
    end;
    if nargout > 4,
      DP22 = minreal ( ss( Ad, B2d, c2, zeros(o2,i2), T ), sqrt(eps), 0);
    end;    
    return
%------- End of SDHINFMOD.M --------- KYuP ----------

%##################################################################    
function [Ad,Bd,Cd,Dd] = ChenFrancis ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 )
%CHENFRANCIS Hinf equivalent model by Chen & Francis.
%
% [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%     Springer-Verlag, Berlin etc., 1995.
%------------------------------------------------------
%	Block sizes
%------------------------------------------------------
    i1 = size(b1,2);
    o1 = size(c1,1);
    n = size(a,1);
%------------------------------------------------------
%	Compute Jinfinity (handle 1-2 block)
%------------------------------------------------------
    av = [a    b2
	      zeros(i2,n+i2)];
    [nv,nv] = size(av);      
	cv = [c1 d12];
	n12 = n + i2;
    ind1 = 1:n12; ind2 = n12+1:2*n12;
	a12 = [    -av'     -cv'*cv
	        zeros(n12)    av   ];
	phi12 = expm ( a12*T );
	Jinf  = - phi12(ind2,ind2)' * phi12(ind1,ind2);
   %------------------------------------------------------
   %	Check using integration of series
   %------------------------------------------------------
   % Jinf1 = matint(av', cv'*cv, av, T, 'single');
   % errRel = norm(Jinf-Jinf1) / norm(Jinf);
   % if errRel > 1.e-6
   %    warning(sprintf('Two procedures for calculating Jinf gave different values: err = %g', errRel ));
   %    keyboard
   % end;
%------------------------------------------------------
%	Form matrices E, X, Y
%------------------------------------------------------
    E = [    -a'   -c1'*c1
         xi2*b1*b1'    a];
    [ne,ne] = size(E); 
    X = [c1 d12]' * [zeros(o1,ne-n) c1];
    Y = [c1 zeros(o1,ne-n)]' * [c1 d12];    
%------------------------------------------------------
%	Compute
%   | P M L |              |-Av'  X  0  |
%   | 0 Q N | = expm ( T * |  0   E  Y  | )
%   | 0 0 R |              |  0   0  Av |
%
%------------------------------------------------------
    aex = [     -av'          X       zeros(nv,nv) 
            zeros(ne,nv)      E            Y
            zeros(nv,nv) zeros(nv,ne)     av ];
    phiEx = expm ( aex*T );
    ind1 = 1:nv; ind2 = nv+1:nv+ne;
    ind3 = nv+ne+1:2*nv+ne;
    M = phiEx ( ind1, ind2 );
    L = phiEx ( ind1, ind3 );
    Q = phiEx ( ind2, ind2 );
    N = phiEx ( ind2, ind3 );
    R = phiEx ( ind3, ind3 );    
%------------------------------------------------------
%	Partition
%       | Q11 Q12 |      | R11 R12 |
%   Q = | Q21 Q22 |, R = |  0   I  | 
%------------------------------------------------------
    Q11 = Q ( 1:n, 1:n );
    Q21 = Q ( n+1:2*n, 1:n );
    R11 = R ( 1:n, 1:n );
    R12 = R ( 1:n, (n+1):nv );
%------------------------------------------------------
%	Compute F = [F1 F2] = [inv(Q11)' 0] M'R
%           Ad = R11 + F1
%           B2d = R12 + F2
%------------------------------------------------------
    Q11inv = inv(Q11);
    n0 = size(M,2) - n;
    F = [Q11inv' zeros(n,n0)] * M' * R;
    Ad = R11 + F(:,1:n);
    B2d = R12 + F(:,(n+1):nv);
%------------------------------------------------------
%	Compute B1d such that B1d*B1d' = gamma^2 Q21 / Q11
%------------------------------------------------------
	BB1d = Q21 / Q11;
    B1d = rectchol(BB1d,size(BB1d,2))'; %sqrtm(BB1d);
    i1d = size(B1d,2);
%------------------------------------------------------
%	Compute [C1d D12d] such that 
%                                   | inv(Q11) 0 |
% [C1d D12d]'*[C1d D12d] = R' * M * |    0     0 | N - R' * L + Jinf
%
%------------------------------------------------------
    RM = R' * M;
    cRM = size(RM,2);
    rN  = size(N,1);  
    U = [ Q11inv  zeros(n,rN-n)
          zeros(cRM-n, rN) ];   
    J = RM * U * N - R' * L + Jinf;    
    o1d = rank(J,1e-10);
	v = rectchol(J,o1d);
	C1d  = v(:, 1:n);
	D12d = v(:, n+1:n12);
%------------------------------------------------------
%	Final discrete Hinf-equivalent model 
%------------------------------------------------------  
    Bd = [ B1d B2d ];
    Cd = [ C1d
           c2 ];
    Dd = [ zeros(o1d,i1d) D12d  
              zeros(o2,i1d+i2) ];          
    return
%------- End of CHENFRANCIS.M --------- KYuP ----------
    
%##################################################################    
function [Ad,Bd,Cd,Dd] = CantoniGlover ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 )
%CANTONIGLOVER Hinf equivalent model by Cantoni and Glover.
%
% [1] M. W. Cantoni, K. Glover, Hinf sampled-data synthesis and 
%     related numerical issues, Automatica, vol. 33, no. 12, 
%     pp. 2233-2241, 1997.  
%------------------------------------------------------
%       K. Yu. Polyakov         06 Feb 2005
%                    		    
%------------------------------------------------------
%   Compute matrix exponential
%------------------------------------------------------
    i1 = size(b1,2);
    o1 = size(c1,1);
    n = size(a,1);
    R = eye(i1); % - d11'*d11/gamma2; % d11 = 0
    S = eye(o1); % - d11*d11'/gamma2; % d11 = 0
    c1d12 = [c1 d12];
    b1zero = [b1; zeros(i2,i1)];
    nbar = n+i2;
    Abar = [a b2;zeros(i2,nbar)];
    Ahat = Abar; % + b1zero*(R\d11')*c1d12/gamma2; % d11 = 0;
    BBhat = b1zero*(R\b1zero')*xi2;
    CChat = c1d12'*(S\c1d12);   
%------------------------------------------------------
%	Form matrices Ehat and Qhat=expm(Ehat*T)
%------------------------------------------------------
    Ehat = [-Ahat' -CChat
             BBhat Ahat ];
    ind1 = 1:nbar;
    ind2 = nbar+1:2*nbar;
%------------------------------------------------------
%       Construct ordered real Schur form
%       and find block sizes
%------------------------------------------------------
    [M,S,vM,n1,n2] = schurdiag ( Ehat, 'i', 1/T );
    if n1 > 0 
       q1 = 1:n1;
       q2 = n1+1:n1+n2;
       M11 = M(ind1,q1); M12 = M(ind1,q2);
       M21 = M(ind2,q1); M22 = M(ind2,q2);
       vM11 = vM(q1,ind1); vM12 = vM(q1,ind2);
       vM21 = vM(q2,ind1); vM22 = vM(q2,ind2);
      
       expS11 = expm(-S(q1,q1)*T);
       expS22 = expm(S(q2,q2)*T);
       ML = [expS11*pinv(M11); null(M11')'];
       IZL = [eye(n1); zeros(nbar-n1,n1)];
       MR = [pinv(vM11')'*expS11 null(vM11)];
       IZR = [eye(n1) zeros(n1,nbar-n1)];
    
       MeM = M12*expS22*vM21;
       OmegaLinv = ML*MeM + IZL*vM11;
       OmegaRinv = M11*IZR + MeM*MR;
       Q11inv = MR / OmegaRinv;
       Q21_Q11inv = (M22*expS22*vM21*MR + M21*IZR) / OmegaRinv;
       Q11inv_Q12 = OmegaLinv \ (ML*M12*expS22*vM22 + IZL*vM12);
    else
%------------------------------------------------------
%	Direct method - without restructuring
%	Partition Qhat accordingly to Ehat
%------------------------------------------------------
       Qhat = expm(Ehat*T);
       Q11 = Qhat(1:nbar,1:nbar);
       Q11inv = inv(Q11);
       Q12 = Qhat(1:nbar,nbar+1:end);
       Q21 = Qhat(nbar+1:end,1:nbar);
       Q22 = Qhat(nbar+1:end,nbar+1:end);
       Q21_Q11inv = Q21 * Q11inv;
       Q11inv_Q12 = Q11inv * Q12;
    end;
%------------------------------------------------------
%	Find matrices of an equivalent discrete-time model
%------------------------------------------------------
    Q11invT = Q11inv';
    Ad = Q11invT(1:n,1:n);
    B2d = Q11invT(1:n,n+1:end);    
    BB1d = Q21_Q11inv(1:n,1:n);
    B1d = rectchol(BB1d,size(BB1d,2))'; %sqrtm(BB1d);    
    CCDD = - Q11inv_Q12;
    C1dD12d = rectchol(CCDD);
    C1d  = C1dD12d(:,1:n);
    D12d = C1dD12d(:,n+1:end);   
%------------------------------------------------------
%	Find an equivalent discrete-time model
%------------------------------------------------------
    i1d = size(B1d,2);
    o1d = size(C1d,1);
    Bd = [B1d B2d];
    Cd = [C1d; c2];
    Dd = [zeros(o1d,n) D12d
           zeros(o2,n+i2)     ];
    return;    
%------- End of CANTONIGLOVER.M --------- KYuP ----------

%##################################################################    
function [Ad,Bd,Cd,Dd] = BamiehPearson ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 )
%BAMIEHPEARSON Hinf equivalent model by Bamieh and Pearson.
%
% [1] Bamieh B.A., Pearson J.B., A general framework for linear periodic systems
%     with applications to $\cal{H}_\infty$ sampled-data control // IEEE
%     Trans. Automat. Contr.}, vol. AC-37, no. 4, pp. 418-435, 1992.
%------------------------------------------------------
%       K. Yu. Polyakov         06 Feb 2005
%                    		    
%------------------------------------------------------
%   Compute matrix exponential
%------------------------------------------------------
    n = size(a,1);
    H = [-a' -xi2*c1'*c1
          b1*b1'  a ];      
%    Gamma = intm ( H, T, 0 ); %expm(H*T);
%    Phi   = intm ( H, T, 1 );
%    Omega = intm ( H, T, 2 );
    H1 = [H eye(2*n); zeros(2*n,4*n)];
    H2 = [H1 eye(4*n); zeros(4*n,8*n)];
    expH2 = expm(H2*T);
    Gamma = expH2(1:2*n,1:2*n);
    Phi   = expH2(1:2*n,2*n+1:4*n);
    Omega = expH2(1:2*n,6*n+1:8*n);   
%------------------------------------------------------
%   Partition matrices
%------------------------------------------------------
    ind1 = 1:n; ind2 = n+1:2*n;
    Gamma11 = Gamma(ind1,ind1); Gamma12 = Gamma(ind1,ind2);
    Gamma21 = Gamma(ind2,ind1); Gamma22 = Gamma(ind2,ind2);
    Phi11 = Phi(ind1,ind1); Phi12 = Phi(ind1,ind2);
    Phi21 = Phi(ind2,ind1); Phi22 = Phi(ind2,ind2);
    Omega11 = Omega(ind1,ind1); Omega12 = Omega(ind1,ind2);
    Omega21 = Omega(ind2,ind1); Omega22 = Omega(ind2,ind2);
%------------------------------------------------------
%   Matrices of the equivalent discrete model
%------------------------------------------------------
    Ad = Gamma22 - Gamma21*(Gamma11\Gamma12);
    BB1d = Gamma21 / Gamma11;
    B1d = rectchol(BB1d,size(BB1d,2))'; %sqrtm(BB1d);
    %Psi = intm ( a, T, 1 );
    B2d = (Phi22 - Gamma21*(Gamma11\Phi12))*b2;    
    CC = - Gamma11 \ Gamma12;
    CD = - Gamma11 \ Phi12 * b2;
    DD = b2'*(Omega12 - Phi11*(Gamma11\Phi12))*b2;
    CCDD = [CC CD;CD' DD];    
    CD = rectchol(CCDD);    
    C1d = CD(:,1:n);
    D12d = CD(:,n+1:end);    
%------------------------------------------------------
%   Finalize the equivalent discrete model
%------------------------------------------------------
    i1d = size(B1d,2);
    o1d = size(C1d,1);
    Bd = [B1d B2d];
    Cd = [C1d; c2];
    Dd = [zeros(o1d,n) D12d
            zeros(o2,n+i2)   ];        
    return;
    
%------- End of BAMIEHPEARSON.M --------- KYuP ----------
    
%#############################################    
function M = intm ( H, t, mult )
%INTM multiple integral of matrix exponential
%
%     M = INTM ( A, T, N )
%
%   Inputs:
%       A - square matrix
%       T - integration limit
%       N - multiplicity of integration
%
%   Outputs:
%	    M = int[0..T] int[0..t(1)] ... int[0..t(N-1)] exp(A*t(N-1)) dt(N-1) ... dt(1) 
%   

% [1] Bamieh B.A., Pearson J.B., A general framework for linear periodic systems
%     with applications to $\cal{H}_\infty$ sampled-data control // IEEE
%     Trans. Automat. Contr.}, vol. AC-37, no. 4, pp. 418-435, 1992.
% Use recursively the formula:
%     int[0..T] exp(A*t) dt = [I 0] * exp([A I;0 0]*t) * [0; I]  
%------------------------------------------------------
%       K. Yu. Polyakov         06 Feb 2005
%                    		    
%------------------------------------------------------
    if mult == 0, M = expm(H*t); return; end;
    n = size(H,1);
    In = eye(n); On = zeros(n);
    L = [In On];
    R = [On; In];
    H1 = [H In; On On];
    M = L * intm ( H1, t, mult-1 ) * R;    
    return    
%------- End of INTM.M --------- KYuP ----------
    
%#######################################################
function [Ad,Bd,Cd,Dd] = MirkinPalmor ( a, b1, b2, c1, c2, d12, T, xi2, i2, o2 )
%MIRKINPALMOR Hinf equivalent model by Mirkin and Palmor
%
% [1] Mirkin L., Tadmor G., Yet another Hinf-discretization //
%     IEEE Trans. Automat. Contr., vol. AC-48, no. 5, pp. 891--894, 2003.
%------------------------------------------------------
%       K. Yu. Polyakov         06 Feb 2005
%                    		    
%------------------------------------------------------
%	Compute matrix exponentials
%------------------------------------------------------
    n = size(a,1);
    Hgamma = [zeros(i2,i2) xi2*d12'*c1   b2'   xi2*d12'*d12
              zeros(n,i2)       a      b1*b1'       b2 
              zeros(n,i2)  -xi2*c1'*c1  -a'   -xi2*c1'*d12    
              zeros(i2,2*i2+2*n) ];
    Hgamma0 = ...
             [zeros(i2,i2) 0*d12'*c1    b2'   0*d12'*d12
              zeros(n,i2)       a      b1*b1'       b2 
              zeros(n,i2)   0*c1'*c1    -a'    0*c1'*d12    
              zeros(i2,2*i2+2*n) ];
    Gamma = expm(Hgamma*T);
    Lam  = expm(Hgamma0*T);
%------------------------------------------------------
%	Partition matrices
%------------------------------------------------------
    ind1 = 1:i2; ind2 = i2+1:i2+n;
    ind3 = i2+n+1:i2+2*n;
    ind4 = i2+2*n+1:2*i2+2*n;
    Lam13 = Lam(ind1,ind3);
    Lam22 = Lam(ind2,ind2);
    Lam23 = Lam(ind2,ind3);
    Lam24 = Lam(ind2, ind4);
    Lam33 = Lam(ind3,ind3);
    Gamma13 = Gamma(ind1, ind3);
    Gamma14 = Gamma(ind1, ind4);
    Gamma22 = Gamma(ind2, ind2);
    Gamma23 = Gamma(ind2, ind3);
    Gamma24 = Gamma(ind2, ind4);
    pinvGamma23 = pinv(Gamma23);
    Gamma32 = Gamma(ind3, ind2);
    Gamma33 = Gamma(ind3, ind3);
    Gamma34 = Gamma(ind3, ind4);
    In = eye(n);    
%------------------------------------------------------
%	Form matrices of discrete-time equivalent model
%------------------------------------------------------
    Ad = Lam22;
    B2d = Lam24;
    B1d = rectchol(Lam23*Lam22', size(Lam22,1))'; %sqrtm(Lam23*Lam22');
    M12 = (In - Lam22'*Gamma33)*pinvGamma23;
    M11 = M12*(Lam22 - Gamma22) - Lam22'*Gamma32;
    M13 = M12*(Lam24 - Gamma24) - Lam22'*Gamma34;
    M22 = Lam33*pinv(Lam23) - Gamma33*pinvGamma23;
    M32 = (Gamma13 - Lam13*Lam22'*Gamma33)*pinvGamma23;
    M33 = M32*(Lam24 - Gamma24) + Gamma14 - Lam13*Lam22'*Gamma34;
   %---------------------------------------------------- 
   %    Symmetrization
   %---------------------------------------------------- 
    %M22 = (M22 + M22') / 2;
    %M33 = (M33 + M33') / 2;
   %---------------------------------------------------- 
   %    Faztorization
   %---------------------------------------------------- 
    CCDD = [M11  M12 M13
            M12' M22 M32'
            M13' M32 M33];
    CD = rectchol(CCDD);
    C1d = CD(:,1:n);
    Dalpha = CD(:,n+1:2*n);
    D12d = CD(:,2*n+1:end);
    D11d = Dalpha*B1d;
%------------------------------------------------------
%	Form an equivalent discrete-time model
%------------------------------------------------------
    i1d = size(B1d,2);
    o1d = size(C1d,1);
    Bd = [B1d B2d];
    Cd = [C1d; c2];
    Dd = [D11d     D12d
          zeros(o2,n+i2)];
      
%------- End of MIRKINPALMOR.M --------- KYuP ----------

%#######################################################
function [Ad,Bd,Cd,Dd] = HayakawaHaraYamamoto ( a, b1, b2, c1, c2, d12, T, xi, i2, o2 )
%HAYAKAWAHARAYAMAMOTO Hinf equivalent model by Hayakawa, Hara and Yamamoto.
%
% [1] Hayakawa Y., Hara S., and Yamamoto Y., Hinf type problem for
%     sampled-data control systems - a solution via minimum energy 
%     characterization // IEEE Trans. Automat. Contr., vol. AC-39, no. 11, 
%     pp. 2278-2284, 1994.
%------------------------------------------------------
%       K. Yu. Polyakov         10 Feb 2005
%                    		    
%------------------------------------------------------
%	Compute matrix exponentials
%------------------------------------------------------
        c1 = c1 * xi;
        d12 = d12 * xi;
        [n,i2] = size(b2);
        Phi = expm( [   a    b2             
                     zeros(i2,n+i2)] * T);
        Ad = Phi(1:n,1:n);
        B2d = Phi(1:n,n+1:end);
        
        Psi = expm([ -a' zeros(n)
                    b1*b1' a] * T);
        ind1 = 1:n; ind2 = n+1:2*n;                       
        W0 = Psi(ind2,ind1) * Psi(ind2,ind2)';
        
        Q = [  -a'  -c1'*c1   -c1'*d12 zeros(n,i2)
             b1*b1'    a         b2    zeros(n,i2)
                     zeros(i2,2*(n+i2))
               b2'  d12'*c1   d12'*d12 zeros(i2)];
        Gamma = expm(Q*T);
        jnd1 = 1:n; jnd2 = n+1:2*n; 
        jnd3 = 2*n+1:2*n+i2; jnd4 = 2*n+i2+1:2*n+2*i2;
        invGamma11 = inv(Gamma(jnd1,jnd1));
        
        W = Gamma(jnd2,ind1) * invGamma11;
        Vcc = - invGamma11 * Gamma(jnd1,jnd2);
        Vcd = - invGamma11 * Gamma(jnd1,jnd3);
        Vdd = Gamma(jnd4,jnd3) + Gamma(jnd4,jnd1)*Vcd;
        M1 = invGamma11 - Phi(ind1,ind1)';
        M2 = Gamma(jnd4,jnd1) * invGamma11 - Phi(ind1,n+1:end)';
        siW = sqrtm(pinv(W));
        N = siW * W0 * siW;
        M = siW * W * siW;
        B1d = sqrtm(W) * sqrtm(N);
        Om = zeros(size(M));
        Omp = zeros(size(M,1), size(Vcd,2));
        R = [Vcc Om   Vcd
             Om   M   Omp
             Vcd' Omp' Vdd];
        R1 = [M1*siW; -sqrtm(N); M2*siW]; 
        CDDd = rectchol(R - R1*R1');     
%------------------------------------------------------
%	Form an equivalent discrete-time model
%------------------------------------------------------
        i1d = size(B1d,2);
        o1d = size(CDDd,1);
        C1d = CDDd(:,1:n);
        D11D12d = CDDd(:,n+1:end);
        Bd = [B1d B2d];
        Cd = [C1d; c2];
        Dd = [D11D12d
             zeros(o2,n+i2)];      

%------- End of HAYAKAWAHARAYAMAMOTO.M --------- KYuP ----------
