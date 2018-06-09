function gamma = d11norm ( A, B, C, D, T, tol )
%D11NORM L2[0,T)-induced norm of the operator D11.
%
%     GAMMA = D11NORM ( A, B, C, D, T, TOL )
%     GAMMA = D11NORM ( SYS, T, TOL )
%
%   Inputs:
%     A, B, C, D - realization of G11
%     SYS - continuous-time LTI-model
%     T   - sampling period
%     TOL - tolerance (default 1e-4)
%
%   Outputs:
%     GAMMA - L2-induced norm of D11
%

% References:
%   [1] Y. Oishi, M. Dahleh, A simple bisection algorithm for the 
%       L2-induced norm of a sampled-data system //  Tech. Rep. 
%       LIDS-P2181, MIT, 1992.
%   [2] T. Chen, B.A. Francis, Optimal Sampled-Data Control Systems, 
%       Berlin etc.: Springer-Verlag, 1995.
%   [3] Modification of the basic idea by K. Polyakov, (C) 2007
%       Gold section search.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Jan-2006 $
%------------------------------------------------------
%       Initialization        
%------------------------------------------------------
        if isa(A,'lti'),
          sys = A;
          T = B;
          if exist('C'), tol = C; end;
        else
          sys = ss(A, B, C, D);
        end;
        [A,B,C,D] = ssdata(minreal(sys));
        n = size(A,1);
        [p,m] = size(sys);
        if ~exist('tol','var'), tol = 1e-4; end;
%------------------------------------------------------
%       Special case - D11 = 0 => gamma11 = 0
%------------------------------------------------------
        [xx,xx,K] = zpkdata(sys);
        if norm(K) < eps 
           gamma11 = 0; 
           return; 
        end;               
%------------------------------------------------------
%       Construct DN matrix using fast sampling
%------------------------------------------------------
        N = 10; 
        % sysN = sdfast ( sys, T, 10, 0, 0 );
        % DN = sysN.d;
        [Af,Bf] = c2d ( A, B, T/N );
        DN = zeros(p*N,m*N);
        Afi = eye(n);
        for i = 0:N-1
          if i == 0
            Dbottom = D;
          else
            Dbottom = [C*Afi*Bf Dbottom]; 
            Afi = Afi * Af; 
          end;
          DN(i*p+1:(i+1)*p,1:(i+1)*m) = Dbottom; 
        end;
        gammaL = norm(DN);
%------------------------------------------------------
%       Calculate gamma using gold section
%------------------------------------------------------
        gammaU = 1.5 * gammaL;
        g = (sqrt(5)-1)/2;
        ind = n+1:2*n;
        while gammaU-gammaL > tol
          gamma1 = g*gammaL + (1-g)*gammaU;
          gamma2 = (1-g)*gammaL + g*gammaU;
          mu1 = muEig ( A, B, C, D, T, gamma1 );
          mu2 = muEig ( A, B, C, D, T, gamma2 );
          if mu1 < mu2, gammaU = gamma2;
          else          gammaL = gamma1; end;
        end;
        gamma = (gammaL + gammaU) / 2;
        
%------- End of D11NORM.M --------- KYuP ----------

%#############################################
function mu = muEig ( A, B, C, D, T, gamma )
%MUEIG Minimal absolute value of svd(Q22)
%
%------------------------------------------------------
%       Check data
%------------------------------------------------------   
        n = size(A, 1);
        p = size(C, 1);
        m = size(B, 2);
%------------------------------------------------------
%       Eliminate D by loop shifting
%------------------------------------------------------   
        R = eye(m) - gamma^2*D'*D;
        S = eye(p) - gamma^2*D*D';
        AA = A + gamma^2*B*(R\D*C);
        CC = C'*(S\C)/gamma^2;
        BB = B*(R\B');
%------------------------------------------------------
%       Compute Mu(gamma)
%------------------------------------------------------   
        %M = [ AA BB
        %     -CC -AA'];
        %ind = n+1:2*n;
        M = [-AA' -CC
              BB  AA];
        ind = 1:n;
        Q = expm(M*T); 
        %mu = min(svd(Q(ind,ind))); 
        mu = min(abs(eig(Q(ind,ind)))); 

%------- End of MUEIG.M --------- KYuP ----------

