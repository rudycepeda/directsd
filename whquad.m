function [K,err,err0,err2] = whquad ( A, B, E, D22, PCancel )
%WHQUAD Wiener-Hopf minimization of frequency-domain quadratic functionals.
%
%     [K,ERR,ERR0,ERR2] = WHQUAD ( A, B, E, D22  )
%     [K,ERR,ERR0,ERR2] = WHQUAD ( A, B, E, D22, PCANCEL  )
%
%   Inputs:
%     A(d), B(d), E(d) - rational functions
%     D22(d) - (discrete) model of plant in the loop
%     PCANCEL - poles that must cancel in L
%
%   Outputs:
%     K - optimal controller as an LTI system
%     ERR  - cost for the optimal system 
%     ERR0 - second term of the cost
%     ERR2 - additional term due to unstable plant
%
%   See also POLQUAD, SSQUAD, POLHINF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Factorization
%------------------------------------------------------        
        T = A.Ts;
        [n,d] = tf2nd ( D22 ); 
        dzpk = sdzpk(d, 1, T);
        Av = minreal(dzpk*(dzpk*dzpk')*dzpk'*A);
        Lam = sfactor(Av);
%------------------------------------------------------
%       Basic controller
%------------------------------------------------------        
        [a0,b0] = dioph ( n, d, poln(1,n.var) );
%------------------------------------------------------
%       Construct L1 and L2
%------------------------------------------------------
        L1 = minreal(dzpk'*dzpk'*sdzpk(B)'/sdzpk(Lam)');
        L2 = minreal(a0*Lam/dzpk);
%------------------------------------------------------
%       Construct L
%------------------------------------------------------
        L = minreal(sumzpk(L1, -L2));
        if norm(L.k) < 1e-10*norm(L1.k), 
           L = zpk([], [], 0, T);
        end;
        if exist('PCancel','var')  &&  ~isempty(PCancel),
          if T ~= 0      
               pCancel = exp(-PCancel*T);
          else pCancel = PCancel; 
          end;
        else
          pCancel = [];  
        end;        
        pL2x = others ( L2.p{1}, L1.p{1} );
        pL = others ( [L1.p{1}; pL2x], pCancel );
        L = setpoles ( L, pL ); 
        [xx,dL] = tf2nd(L);
        [dLp,xx,dL0] = factor ( dL );
        dLp = dLp*dL0;        
%------------------------------------------------------
%       Construct Lplus and Psi
%------------------------------------------------------
        [pL,Lsp] = improper ( L );
        [Lp,Lm] = separss ( Lsp, dLp.z );
        Lp = Lp + tf(pL,1,T);
          %[nMinus,nPlus] = dioph ( dLp, dLm, nL );
          %Lp = zpk ( nPlus, dLp, T );
          %Lm = zpk ( nMinus, dLm, T );
        Psi = minreal(Lp/Lam);
%------------------------------------------------------
%       Construct controller
%------------------------------------------------------
        Kn = zpk(poln(a0,'z'),1,T) + dzpk*Psi;
        Kd = zpk(poln(b0,'z'),1,T) - zpk(n,1,T)*Psi;
        K = minreal(Kn/Kd);
        [nK,dK] = tfdata ( K, 'v' );
        nK = striplz(nK);
        dK = striplz(dK);
        K = zpk(tf(nK,dK,K.Ts));
          %a0,b0,Psi
          %Kn = tf(Kn)
          %Kd = tf(Kd)
        if T ~= 0, K = K'; end;
%------------------------------------------------------
%       Unstable part of L
%------------------------------------------------------
        if T == 0, Lm = Lm'; end;
        errLm = norm(Lm)^2;
%------------------------------------------------------
%       Compute int_Gamma (E - B*B'/A) dz/z
%------------------------------------------------------
        BBA = minreal(B*B'/A);
         if T == 0
           X = sumzpk(E, -BBA);
           F = sfactor ( X );
           err0 = norm(F)^2;
        else
           err0 = intz2b (E, -BBA); 
        end;
        err = errLm + err0;
%------------------------------------------------------
%       Additional term due to unstable plants
%------------------------------------------------------
        if nargout > 3
          [xx,L2] = improper(L2);
          [xx,dL2] = tf2nd(L2);
          [xx,dL2m,xx] = factor ( dL2 );
          if deg(dL2m) > 0
            L2m = separss(L2, dL2m.z, 'infu');            
            err2 = norm(L2m)^2;
          else
            err2 = 0;  
          end;
        end;

%------- End of WHQUAD.M --------- KYuP ----------           
        
