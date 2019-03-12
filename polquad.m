function [K,err,err0,err2] = polquad ( A, B, E, D22, PCancel, Sigma )
%POLQUAD Polynomial minimizition of a quadratic functional.
%
%     [K,ERR,ERR0,ERR2] = POLQUAD ( A, B, E, D22  )
%     [K,ERR,ERR0,ERR2] = POLQUAD ( A, B, E, D22, PCANCEL  )
%
%   Inputs:
%     A(d), B(d), E(d) - rational functions in 's' or 'zeta'
%     D22(d) - (discrete) model of plant in the loop
%     PCANCEL - poles that must cancel in L
%
%   Outputs:
%     K - optimal (discrete) controller as an LTI system
%     ERR  - cost for the optimal system 
%     ERR0 - second term of the cost
%     ERR2 - additional term due to unstable plant
%
%   See also SSQUAD, WHQUAD, QUADERR, POLHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Take account of Sigma
%------------------------------------------------------        
        if ~exist('PCancel','var')
            PCancel = []; 
        end
        if exist('Sigma','var')
           %A = minreal(A*Sigma);
           %B = minreal(B*Sigma);
           %E = minreal(E*Sigma);
          [nSigma,dSigma] = tf2nd(sfactor(Sigma));   
        else
          nSigma = 1;   
          dSigma = 1;   
        end
%------------------------------------------------------
%       Construct YN and YD
%------------------------------------------------------        
        T = A.Ts;
        [n,~,~,~,Lam,na,da,nb,db,dLp,dLm,dL0,nL,YN,YD,~,~,~,L2] = ...
             ynyd ( A, B, D22, PCancel );
%------------------------------------------------------
%       Solve polynomial equation
%------------------------------------------------------
       %-----------------------------------------------
       %       Open-loop system
       %-----------------------------------------------
        if norm(n) < eps % 
          AN = dLp; BN = dLm; CN = nL;
         %-----------------------------------------------
         %   Extract common divisors of AN, BN and CN
         %-----------------------------------------------
          %[AN,BN,CN] = triple ( AN, BN, CN );
          [AN,BN,gABN] = coprime(AN, BN);
          [CN,remN] = mrdivide(CN, gABN);
          if norm(remN) > 1e-3*norm(CN)
             warning('Incorrect division of CN: err=%g',norm(remN)); 
          end
         %-----------------------------------------------
         %   Solve polynomial equation
         %-----------------------------------------------
          [AN,CN,gN] = coprime(AN, CN);
          [Pi,N] = dioph(AN, BN, CN);        
          K = minreal(sdzpk(gN*N, dLp, A.Ts)/sdzpk(Lam), 1e-3);
          K = sdzpk(K);
       %-----------------------------------------------
       %       Closed-loop system
       %-----------------------------------------------
        else
          AN = dLp*na*dSigma;   BN = dLm;  CN = nSigma*YN;
          AD = - dLp*nb*dSigma; BD = dLm;  CD = nSigma*YD;
         %-----------------------------------------------
         %   Extract common divisors of AN, BN and CN; AD, BD and CD
         %-----------------------------------------------
          %[AN,BN,CN] = triple ( AN, BN, CN );
          %[AD,BD,CD] = triple ( AD, BD, CD );
          [AN,BN,gABN] = coprime(AN, BN);
          [CN,remN] = mrdivide(CN, gABN);
          if norm(remN) > 1e-3*norm(CN)
             warning('Incorrect division of CN: err=%g',norm(remN)); 
          end
          [AD,BD,gABD] = coprime(AD, BD);
          [CD,remD] = mrdivide(CD, gABD);
          if norm(remD) > 1e-6*norm(CD)
             warning('Incorrect division of CD: err=%g',norm(remD)); 
          end
         %-----------------------------------------------
         %   Extract possible divisor of CD at the unit circle
         %-----------------------------------------------
          gD0 = 1; 
          if deg(dL0) > 0
            [CDx,err] = mrdivide(CD, dL0);
            if norm(err) < 1e-3*norm(CD)
               AD = AD / dL0;
               CD = CDx;
               gD0 = dL0;
            end
          end
         %-----------------------------------------------
         %   Solve polynomial equations
         %-----------------------------------------------
          [AN,CN,gN] = coprime(AN, CN);
          [AD,CD,gD] = coprime(AD, CD);
          [Pi,N,D] = diophsys(AN, BN, CN, AD, BD, CD);
          if norm(D) < eps*norm(N)
             K = Inf; err = Inf; err0 = 0; err2 = 0;  
             return; 
          end
          K = minreal(sdzpk(db*gN*N, da*gD0*gD*D, A.Ts), 1e-3);
          K = sdzpk(K);
        end
        if T ~= 0
            K = K'; 
        end
%------------------------------------------------------
%       Return if only controller needed
%------------------------------------------------------
        if nargout < 2
            return
        end
%------------------------------------------------------
%       Unstable part of L
%------------------------------------------------------
        Lm = sdzpk(Pi, dLm, A.Ts);
        if T == 0
            Lm = Lm';
        end
        errLm = norm(Lm)^2;
%------------------------------------------------------
%       Compute int_Gamma (E - B*B'/A) dz/z
%------------------------------------------------------
        BBA = minreal((sdzpk(B)*sdzpk(B)')/A);
        if T == 0
           X = sumzpk(E, -BBA);
           F = sfactor ( X );
           err0 = norm(F)^2;
        else
           err0 = intz2b (E, -BBA); 
        end
        err = errLm + err0;
%------------------------------------------------------
%       Additional term due to unstable plants
%------------------------------------------------------
        if nargout > 3
          [~,L2] = improper(L2);
          [~,dL2] = tf2nd(L2);
          [~,dL2m,~] = factor ( dL2 );
          if deg(dL2m) > 0
            L2m = separss(L2, dL2m.z, 'infu');            
            err2 = norm(L2m)^2;
          else
            err2 = 0;  
          end
        end

%------- End of POLQUAD.M --------- KYuP ----------           
        
