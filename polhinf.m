function [K,lam,Sigma] = polhinf ( A, B, E, D22, PCancel )
%POLHINF Polynomial solution to an Hinf-problem.
%
%     [K,LAM,SIGMA] = POLHINF ( A, B, E, D22  )
%     [K,LAM,SIGMA] = POLHINF ( A, B, E, D22, PCANCEL )
%
%   Inputs:
%     A(d), B(d), E(d) - discrete-time functions
%     D22(d) - discrete model of plant in the loop
%     PCANCEL - poles that must cancel in L
%
%   Outputs:
%     K   - Hinf-optimal discrete controller as an LTI system
%     LAM - cost function for the optimal system
%     SIGMA - optimal weighting function for the optimal system
%
%   See also DHINF, SDAHINF, SDTRHINF, POLQUAD.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Construct YN and YD
%------------------------------------------------------        
        if ~exist('PCancel','var'), PCancel = []; end;
        T = A.Ts;
        [n,d,a0,b0,Lam,na,da,nb,db,dLp,dLm,dL0,nL,YN,YD,...
         dLCancel,L] = ...
             ynyd ( A, B, D22, PCancel );
%------------------------------------------------------
%       Form Z
%------------------------------------------------------
        [Z,phi2,dz] = zterm ( A, B, E, dLCancel );
        dz2 = dz * dz';
%------------------------------------------------------
%       Construct g-factor
%------------------------------------------------------
        if 0 
           V2 = minreal(dL*Lam/dzpk);
           nV2 = tf2nd ( V2 );
           zV2s = extrpair ( nV2.z, 'd' );
           [xx,g] = others(zV2s, roots(dz));
        else
           g = []; % better numerically!
        end;
        g = poln(g, ['r' d.var]);
%------------------------------------------------------
%       Construct \eta and \chi
%------------------------------------------------------
        Y = zpk(dz, delzero(dLm), Z.Ts);
        FY = minreal( sfactor(Y*Y'), 1e-4 );
        [chi,eta] = tf2nd ( FY );
%------------------------------------------------------
%       Extract common factors
%------------------------------------------------------
        b0 = dLm; q0 = eta*dLp; c0 = chi*nL; r0 = g;
        [b0,q0,c0] = triple ( b0, q0, c0 );
        [r0,q0,c0] = triple ( r0, q0, c0 );
        [q0,c0,gQ] = coprime ( q0, c0 );
%------------------------------------------------------
%       Find optimal lam using single equation
%------------------------------------------------------
        ws = warning; warning off;
        [lam,P,Q,sigma,count,err0] = hinfred ( dz2, phi2, b0, q0, c0, r0 );
        if isinf(lam)
           disp('HINFRED failed'); 
           [lam,P,Q,sigma,xx,xx,err0] = hinffiter  ( dz2, phi2, b0, q0, c0, r0 );
           if isinf(lam)
             disp('HINFFITER failed'); 
             [lam,P,Q,sigma,count,err0] = hinfbisec  ( dz2, phi2, b0, q0, c0, r0 );
             if isinf(lam), 
                disp('HINFBISEC failed'); 
                error('Numerically unstable AHinf-problem!'); 
             end;
           end;
        end;
        warning(ws);            
%        [lam,P,Q,sigma,count,err0] = hinfnehari ( dz2, phi2, b0, q0, c0 );
%        lam, P, sigma, err0
%------------------------------------------------------
%       Check non-generic case
%------------------------------------------------------
        nonGen = iscell(P);
        if nonGen
           P1 = P{2}; P = P{1};
        end;
%------------------------------------------------------
%       Check solution of reduced polynomial equation
%------------------------------------------------------
        tol = 1e-3;        
        p10 = b0*r0*Q + P*sigma*q0;
        p20 = recip(P)*c0;
        errP120 = norm(p10-p20);
        if errP120 > tol*norm(p10)
           disp(sprintf('Incorrect solution: errP120 = %g',errP120)); 
        end;
%------------------------------------------------------
%       Check solution of initial polynomial equation
%------------------------------------------------------
        Q = Q * gQ;        
        p1 = dLm*Q + P*sigma*eta*dLp;
        p2 = recip(P)*chi*nL;
        errP12 = norm(p1-p2);
        if errP12 > tol*norm(p1)
           disp(sprintf('Incorrect solution: errP12 = %g',errP12)); 
        end;
%------------------------------------------------------
%       Form optimal weight Sigma
%------------------------------------------------------
        if 1  % nargout > 2,
          Fsigma = zpk(recip(P)*chi, sigma*eta, T);
          Sigma = Fsigma * Fsigma';
        end;
%------------------------------------------------------
%       Check equality L*FSi = W + Y
%------------------------------------------------------       
if 1
        lamZ = sumzpk(-Z,lam^2);
        LFsigma = minreal(L*Fsigma,tol);
        W = zpk(Q, dLp*sigma*eta,T);
        Y = zpk(P, dLm, T);
        LFsigma1 = sumzpk(W,Y);
        dLFsigma = sumzpk(LFsigma,-LFsigma1);
        errLFSi = norm(tfdata(dLFsigma,'v'))/norm(tfdata(LFsigma,'v')); 
                  %max(bode(dLFsigma));
        if errLFSi > tol,
           disp(sprintf('Incorrect solution: errLFSi = %g',errLFSi)); 
        end;
end;        
%------------------------------------------------------
%       Check equality YY = FSi*FSi'(lam^2 - Z)
%------------------------------------------------------       
if 1
        YY = minreal(Y*Y');
        ZZ = minreal(Sigma*lamZ);
        YYmZZ = sumzpk(YY, -ZZ);
        errYZ = max(bode(YYmZZ));
        if errYZ > tol,
           disp(sprintf('Incorrect solution: errYZ = %g',errYZ)); 
        end;
end;        
%------------------------------------------------------
%       Check equality X = lam^2
%------------------------------------------------------              
if 1
        PsiLam = minreal(sdzpk(Q,dLp*recip(P)*chi,T), tol);
        W  = minreal(W,tol);
        W1 = minreal(PsiLam*Fsigma,tol);
        dW = sumzpk(W,-W1);
        errW = norm(tfdata(dW,'v'))/norm(tfdata(W,'v')); %max(bode(dW));
        if errW > tol,
           disp(sprintf('Incorrect solution: errW = %g',errW)); 
        end;        
        U = minreal(Y/Fsigma, 1e-4);
        U = sdzpk(U);
        X = sumzpk(minreal(U*U'),Z);
        dX = sumzpk(X,-lam^2);
        errX = max(bode(dX))/max(bode(X));
        if errX > 1e-3*lam,
           disp(sprintf('Incorrect solution: errX = %g',errX)); 
        end;        
end;        
%------------------------------------------------------
%       Find N and D by division
%------------------------------------------------------
        U = sigma*eta*dLp;
        AN =   U*na; BN = - chi*YN;
        AD = - U*nb; BD = - chi*YD;
        CN = dLm*g; CD = CN;                

        [AN,BN,CN] = triple ( AN, BN, CN );
        [AD,BD,CD] = triple ( AD, BD, CD );
        [AN,BN,gN] = coprime(AN, BN);
        [AD,BD,gD] = coprime(AD, BD);

        N0 = - recip(P)*BN - P*AN;
        [N,errN] = mrdivide(N0, CN);
        if norm(errN) > 1e-6*norm(N0),
           warning('Incorrect division: errN = %g',norm(errN)); 
        end;

        D0 = - recip(P)*BD - P*AD; 
        [D,errD] = mrdivide(D0, CD);
        if norm(errD) > 1e-6*norm(N0),
           warning('Incorrect division: errD = %g',norm(errD)); 
        end;
%------------------------------------------------------
%       Form controller
%------------------------------------------------------
        K = minreal(sdzpk(db*gN*N, da*gD*D, T), 1e-3);
        K = sdzpk(K);
        K = minreal(K', 1e-3);        
%------------------------------------------------------
%       Form second non-generic controller
%------------------------------------------------------
        if nonGen
          AN = - AN; 
          N0 = - recip(P1)*BN - P1*AN;
          [N,errN] = mrdivide(N0,CN);
          if norm(errN) > 1e-6*norm(N0),
             warning('Incorrect division: errN = %g',norm(errN)); 
          end;
          AD = - AD; 
          D0 = - recip(P1)*BD - P1*AD; 
          [D,errD] = mrdivide(D0,CD);
          if norm(errD) > 1e-6*norm(N0),
             warning('Incorrect division: errD = %g',norm(errD)); 
          end;
          %Edited by RCG
          nn = db*gN*N;
          dd = da*gD*D;
          K1 = sdzpk(nn,dd,T);
          %K1 = minreal(K1);
          K1 = minreal(K1', 1e-3); 
          K = {K K1};
        end;
        
%------- End of POLHINF.M --------- KYuP ----------           

       
