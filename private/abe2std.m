function [dsys,X] = abe2std ( A0, A1, B, E, G22 )
%ABE2STD State-space model for a freqency-domain functional.
%
%     [SYS,X] = ABE2STD ( A0, A1, B, E, D22 )
%
%   Inputs:
%     A0(d), A1(d), B(d), E(d) - discrete-time functions
%     D22(d) - discrete model of plant in the loop
%
%   Outputs:
%     SYS - discrete-time state-space standar
%     X   - additional term, X = E - G11'*G11  
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if any(size(G22) ~= 1)
           error('Function implements only scalar controller design'); 
        end;
%------------------------------------------------------
%       Factorizations
%------------------------------------------------------        
        if 0 
          A0v = minreal(dzpk*dzpk'*A0);
          A1v = minreal(dzpk*dzpk'*A1);
        else
          [nLam0,dLam0] = tf2nd(sfactor(A0, 'd'));
          Lam0 = zpk(recip(nLam0), dLam0, A0.Ts);       
          [nLam1,dLam1] = tf2nd(sfactor(A1, 'd'));
          Lam1 = zpk(recip(nLam1), dLam1, A1.Ts);
        end;
%------------------------------------------------------
%       G12 and G21
%------------------------------------------------------        
        if 0
          G12 = minreal(Lam0 / dzpk);
          G21 = minreal(Lam1 / dzpk);
        else
          G12 = Lam0; 
          G21 = Lam1; 
        end;
%------------------------------------------------------
%       G11
%------------------------------------------------------
        if 0 
           G110 = - minreal(dzpk'*dzpk'*B'/Lam0'/Lam0');
           [n11,d11] = tf2nd ( G110 );
           p11d = roots(d11);
           p22d = G22.p{1};
           p0 = others ( p11d, p22d );
           pDelete = p0(abs(p0) < 1);
           [Gu,G11] = separss ( G110, pDelete, 'infu' );
        else
           G11 = - minreal(B/G12/G21);
           [p11,G11] = improper ( G11, 'p' );
           G11 = z2zeta(G11);
        end;
%------------------------------------------------------
%       Full standard system
%------------------------------------------------------
        dsys = z2zeta(zpk([G11 G12; G21 -G22]));
        dsys = minreal(ss(dsys), sqrt(eps), 0);
        X = minreals(sumzpk(E,-G11'*G11), 1e-3);

%------- End of ABE2STD.M --------- KYuP ----------           
        
