function [n,d,a0,b0,Lam,na,da,nb,db,dLp,dLm,dL0,nL,YN,YD,...
          dLCancel,L,L1,L2] = ...
         ynyd ( A, B, D22, PCancel )
%YNYD Find auxiliary polynomials for polynomial sampled-data design.
%
%     [n,d,a0,b0,Lam,na,da,nb,db,dLp,dLm,dL0,nL,YN,YD,dLCancel] = ...
%         YNYD ( A, B, D22, PCancel )
%
%   Inputs:
%     A(d), B(d) - rational functions in 's' or 'zeta'
%     D22(d)  - (discrete) model of plant in the loop
%     PCANCEL - poles that must cancel in L
%
%   Outputs:
%     (n,d)   - polynomials forming D22
%     (a0,b0) - basic controller
%     Lam     - spectral factor of A
%     (na,da), (nb,db) - pairs of coprime polynomials
%      dLp, dLm, dL0 - stable, unstable and neutral factors of DL
%      nL - numerator of L
%      YN, YD - polynomials in the right parts of equations
%      dLCancel - cancelled polynomial in L
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 29-Sep-2006 $
%------------------------------------------------------
%       Factorization
%------------------------------------------------------        
        T = A.Ts;
        [n,d] = tf2nd ( minreal(zpk(D22)) ); 
        dzpk = zpk(d, 1, T);
        Av = dzpk*(dzpk*dzpk')*dzpk'*A;
        Lam = minreal( sfactor(Av) );
%------------------------------------------------------
%       Basic controller
%------------------------------------------------------        
        [a0,b0] = dioph ( n, d, poln(1,n.var));
%------------------------------------------------------
%       Construct L1 and L2
%------------------------------------------------------
        L1 = minreal(dzpk'*dzpk'*B'/Lam');
        if norm(a0) > eps 
             L2 = minreal(a0*Lam/dzpk);        
        else L2 = zpk([], [], 0, T); 
        end;        
%------------------------------------------------------
%       Construct L
%------------------------------------------------------
        L = minreal(sumzpk(L1, -L2));
        if norm(L.k) < 1e-10*norm(L1.k), 
          L = zpk([], [], 0, T);
        end;
%------------------------------------------------------
%       Correct poles of L
%------------------------------------------------------
        pL2x = others ( L2.p{1}, L1.p{1} );
        pL = [L1.p{1}; pL2x];
        if exist('PCancel','var') && ~isempty(PCancel),
          if T ~= 0      
               pCancel = exp(-PCancel*T);
          else pCancel = PCancel; 
          end;
          pL = others ( pL, pCancel );
          dLCancel = poln(pCancel, 'rz');
        else
          dLCancel = poln(1, 'z');
        end;       
        L = setpoles ( L, pL );        
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
           % L = minreal(L, 1e-3); % May be incorrect!!!
        [nL,dL] = tf2nd(L);
        [dLp,dLm,dL0] = factor ( dL );
        dLp = dLp*dL0;
        
        dLLam = minreal(zpk(dLp,1,T)*Lam);
        [na,da] = tf2nd(minreal(dzpk/dLLam));
        [nb,db] = tf2nd(minreal(zpk(n,1,T)/dLLam));        
%------------------------------------------------------
%       Construct YN and YD
%------------------------------------------------------
        R1 = minreal(zpk(dL,1,T)*L1);
        YN = minreal(zpk(na,1,T)*R1);
        [nYN,dYN]  = tf2nd ( YN );        
        [YN,remN] = mrdivide ( nYN, dYN );
        if norm(remN) > 1e-3*norm(nYN)
           warning(sprintf('Incorrect division of YN: err=%g',norm(remN))); 
        end;                
        YD1 = minreal(zpk(dLm*db,1,T)/dzpk);
        YD2 = minreal(R1*zpk(nb,1,T));
        YD  = minreal(sumzpk(YD1,-YD2),1e-4);        
        [nYD,dYD] = tf2nd ( YD );
        if norm(nYD) < 1e-3*(norm(tf2nd(YD1))+norm(tf2nd(YD2)))
          nYD = nYD*0;  
        end;
        [YD,remD] = mrdivide ( nYD, dYD );
        if norm(remD) > 1e-3*norm(nYD)
           warning(sprintf('Incorrect division of YD: err=%g',norm(remD))); 
        end;        
        
%------- End of YNYD.M --------- KYuP ----------           
