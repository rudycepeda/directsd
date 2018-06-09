function [r,k,rts] = polfactz ( p, tol )
%POLFACTZ Factorize a symmetric quasipolynomial.
%         Newton iteration method by Z. Vostry.
%
%     [F,K,R] = POLFACTZ ( P, TOL )
%
%   Inputs:
%     P - symmetric quasipolynomial such that P(z)=P(1/z) 
%     TOL  - tolerance (default tol = 1e-6)
%
%   Outputs:
%     F - stable factor (with roots OUTSIDE the unit circle)
%         such that F*F~ = P
%     K - leading coefficient of F
%     R - roots of F (K-R representation is more reliable for
%         polynomials with roots on the unit circle)
%

% [1] Z. Vostry, New Algorithm for Polynomial Spectral Factorization 
%     with Quadratic Convergence I, Kybernetika, 11, 6, 415-422, 1975.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1e-6; end;
        if norm(p - fliplr(p)) > tol, 
          error('Quasipolynomial must be symmetric'); 
        end;  
%------------------------------------------------------
%       Extract roots at z=1 and z=-1
%------------------------------------------------------
        r0 = 1;
        rts = [];
        for i=1:2
           if i == 1, q = 1; else q = -1; end;
           d = conv([1 -q], [-q 1]);
           while 1
              [p0,rem] = poldiv(p, d);    
              if abs(rem) > tol, break; end;
              p = p0;
              r0 = conv(r0, [1 -q]);
              rts = [rts; q];              
           end;
        end;    
%------------------------------------------------------
%       Newton iteration loop
%------------------------------------------------------
        n = ceil(length(p) / 2);
        B = 2*p(1:n).';        
        r = fliplr(eye(1,n));
        while 1
           A = toep(fliplr(r),n,n) + hank(r,n,n); 
           vecS = A \ B;
           s = fliplr(vecS.');
           rOld = r;
           r = (r + s) / 2;
           if abs(rOld - r) < tol, break; end;
        end;    
%------------------------------------------------------
%       Final results
%------------------------------------------------------
        r = striplz(r);
        k = r(1);
        rts = [rts; roots(r)];
        r = conv(r, r0);

%------- End of POLFACTZ.M --------- KYuP ----------
%#################################################    
function [quo,rem] = poldiv ( p, d )
%-----------------------------------------------------
%
% POLDIV more reliable polynomial division instead deconv
%
%-----------------------------------------------------
    lenP = length(p);
    lenD = length(d);
    if lenP < lenD,
       quo = 0;
       rem = p;
       return;
    end;    
    lenQ = lenP - lenD + 1;
    lenR = lenD - 1;    
    A = [toep(d,lenP,lenQ) eye(lenP,lenR)];
    B = toep(p, lenP, 1);   
    X = A \ B;
    quo = fliplr(X(1:lenQ).');
    rem = fliplr(X(lenQ+1:end).');
    
%------- End of DIVPOL.M --------- KYuP ----------
    
