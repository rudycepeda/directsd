function [Q,E,Emin] = polopth2 ( V, W, n )
%POLOPTH2 Solution to polynomial H2-minimization problem. 
%          ||V(d) -  W(d)Q(d)||_{H2} -> min
%          
%    [Q,E,EMIN] = POLOPTH2 ( V, W, n )
%                                       
%  Inputs:
%	 V, W - stable LTI discrete models in d=z^{-1}
%    n    - degree of desired Q(d) = q_n d^n + ... + q_1 d + q_0
%
%  Outputs:
%    Q - optimal polynomial solution or empty if none exists
%    E - minimal attainable H2-norm for deg(Q)=n
%    EMIN - minimal attainable H2-norm for deg(Q)->Inf
%
%   See also MODSDH2, MODSDL2, POLQUAD, WHQUAD.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Construct polynomials for W, V and W^*
%------------------------------------------------------
        Wz = W';
        [nW,dW]   = tfnumden(W);
        [nV,dV]   = tfnumden(V);
        [nWz,dWz] = tfnumden(Wz);        
%------------------------------------------------------
%       Construct polynomials for polynomial equation
%------------------------------------------------------
        pQ = conv(nWz,conv(nW,dV));
        pP = - conv(dW,dV);
        pT = - [dWz zeros(1,n+1)];
        c = conv(nWz, conv(nV,dW));        
%------------------------------------------------------
%       Determine degrees of unknown polynomials
%------------------------------------------------------
        %degP = deg(pT) - 2 - n;
        degP = deg(dWz) - 1;
        degQ = n;
        maxDeg = max([deg(pQ)+degQ deg(pP)+degP deg(c)]);
        degT = maxDeg - deg(pT);        
%------------------------------------------------------
%       Form linear system of equations
%------------------------------------------------------
        AQ = toep ( pQ, maxDeg+1, degQ+1 );
        AP = toep ( pP, maxDeg+1, degP+1 );
        AT = toep ( pT, maxDeg+1, degT+1 );
        A = [AQ AP AT];
        B  = toep ( c, maxDeg+1, 1 );       
%------------------------------------------------------
%       Solve linear system of equations
%------------------------------------------------------
        [flag,X] = linsys ( A, B );
        if flag == 0
          ws = warning; warning off; 
          X = A \ B;   % it cannot be worse
          warning(ws);
        end;
        indQ = 1:degQ+1;
%        indP = degQ+2:degQ+1+degP+1;
%        indT = degQ+1+degP+2:degQ+1+degP+1+degT+1;
%------------------------------------------------------
%       Extract polynomials
%------------------------------------------------------
        Q = flipud(X(indQ)).';
%        P = flipud(X(indP)).';
%        T = flipud(X(indT)).';       
%        wVW = dimpresp ( V/W, 10, 'd' )
%        R = minreal(tf(P, nWz, W.Ts) / W);
%        wR = dimpresp ( R, 10, 'd' )
%        wX = wVW + wR
%        Q        
%------------------------------------------------------
%       Limiting case deg(Q)->Inf
%------------------------------------------------------
        if nargout > 2
           rtsNw = roots(nW);
           ind = find(abs(rtsNw) <= 1);
           if isempty(ind)
              Emin = 0;
           else
              mu = sfactz ( nWz, 'zeta' );
              pi0 = dioph ( dV, - poly(rtsNw(ind)), - conv(nV,mu)); 
              R = tf ( pi0, mu, W.Ts );
              Emin = norm(R')^2;
%              nWzm = deconv( nWz, mu );              
%              piInf = conv( nWzm, pi0 );
           end;     
        end;
%------------------------------------------------------
%       Final results
%------------------------------------------------------
        zpkQ = zpk(tf(Q, 1, V.Ts));
        E = norm(minreal(V - W*zpkQ)')^2;
%        Q = poln(Q, 'z');

%------- End of POLOPTH2.M --------- KYuP ----------
%##################################################
function [num,den] = tfnumden ( F )
%
% TFNUMDEN find numerator and denominator of a scalar transfer function
%          and strip leading zeros
%
%------------------------------------------------------
%       K. Yu. Polyakov         12 Dec 2004
%                                          
%------------------------------------------------------
        F = tf(F);
        if max(size(F)) > 1, error('Operation is applicable to scalar functions only'); end;
        [num,den] = tfdata ( F, 'v' );
        num = striplz ( num );
        den = striplz ( den );
        
%------- End of TFNUMDEN.M --------- KYuP ----------
