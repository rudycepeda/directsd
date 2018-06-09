function [Z,phi2,dz] = zterm ( A, B, E, dLCancel ) 
%ZTERM Calculate Z = E - B*B'/A with cancellations.
%
%     [Z,NZ2,DZ] = ZTERM ( A, B, E, DL0 )
%
%   Inputs:
%     A, B, E - matrix coefficients in quadratic functionals
%     DL0 - polynomial containing all cancelled poles of L
%
%   Outputs:
%     Z   - cancelled function Z = E - B*B'/A
%     NZ2 - numerator of Z, symmetric quasipolynomial
%     DZ  - stable factor of denominator Z
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Separate cancellable part of E
%------------------------------------------------------
     T = E.Ts;
     [~,dE] = tf2nd ( E );
     [~,~,dCancel] = coprime(dE, dLCancel*dLCancel');        
%------------------------------------------------------
%       Direct way without cancellations
%------------------------------------------------------
     BBA = minreals(sdzpk(B*B'/A));
     Z = sumzpk(E,-BBA);
     if norm(Z.k) < eps  ||  norm(Z.k) < 1e-8*norm(E.k)
        Z = zpk([], [], 0, Z.Ts);  
     end;
     Z = symmetr(Z);
     Z = minreals(sdzpk(Z), 1e-3);
     [~,zCancel] = others(Z.p{1}, dCancel.z);
%------------------------------------------------------
%       Special construction with cancellation
%------------------------------------------------------
     if ~isempty(zCancel)
        [E0,E1] = separss ( E, dCancel.z, 'infu' );                
       %------------------------------------------------------
       %       Separate cancellable part of dA/(dB*dB')
       %------------------------------------------------------        
        [nB,dB] = tf2nd(B);
        [nA,dA] = tf2nd(A);
        [nD,dD] = tf2nd(minreal(zpk(dA,dB*dB',T)));
        [nD,znD] = delzero(nD);
        [dD,zdD] = delzero(dD);
        DD = zpk(nD, dD, T);
        [DD0,DD1] = separss ( DD, dCancel.z, 'infu' );             
        [nD0,dD0] = tf2nd ( DD0 ); 
       %------------------------------------------------------
       %       Separate cancellable part of B*B'/A
       %------------------------------------------------------        
        [nN,dN] = tf2nd(minreal(zpk(nB*nB',nA,T)));
        [nN,znN] = delzero(nN);
        [dN,zdN] = delzero(dN);
        NN = zpk(nN, dN, T);
        zerosND = znD + znN  - zdD - zdN;
        if zerosND > 0, nN = nN * z^zerosND; end;
        if zerosND < 0, dN = dN * z^(-zerosND); end;
        [m0,m1] = dioph(dN, dD0, nN*nD0);
        X0 = zpk(m0, dD0, T);
        X1 = zpk(m1, dN, T);
       %------------------------------------------------------
       %       Construct Z = E-B*B'/A
       %------------------------------------------------------
        diff = minreal(sumzpk(E0, -X0));
        if norm(diff.k) < sqrt(eps), diff = 0; end;
        Z0 = minreal(sumzpk(X1, NN*DD1), 1e-3);
        Z0 = sumzpk(Z0, Z0')/2;
        Z = sumzpk(E1, -Z0);
     end;
%------------------------------------------------------
%    Extract poles at z=1, which should be cancelled
%------------------------------------------------------
     tol = sqrt(eps);
     [zZ,pZ,kZ] = zpkdata(Z, 'v');
     ind = find(abs(pZ-1) < tol);
     if ~isempty(ind)
       for i=1:length(ind)
         try pZ = remove ( pZ, 1 );  catch end;  
         try zZ = remove ( zZ, 1 );  catch end;
       end;     
       Z.z{1} = zZ; Z.p{1} = pZ;
     end;
%------------------------------------------------------
%       Find dz and phi
%------------------------------------------------------
     [nZ,dZ] = tf2nd(Z);
     [dZ,zerosD] = delzero(dZ);
     ell = floor(deg(dZ)/2);
     mid = zerosD + ell;
     cut0 = 2*mid - deg(nZ);
     if cut0 > 0, 
        nZc = nZ.coef;
        nZc(end-cut0+1:end) = 0;
        nZ = poln(nZc,'z');
     end;
     try   K = 1;  dz = sfactor(dZ); 
     catch K = -1; dz = sfactor(-dZ); end;
%------------------------------------------------------
%    Symmetrize phi2
%------------------------------------------------------
     phi2 = delzero(K*nZ);        
     phi2 = phi2*z^(-deg(phi2)/2);
%------------------------------------------------------
%    Make dZ monic
%------------------------------------------------------
     phi2 = phi2 / dz.coef(1)^2;
     dz = dz / dz.coef(1);
        
%------- End of ZTERM.M --------- KYuP ----------
