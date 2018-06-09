function [lam,P,Q,sigma,count,err] = hinfbisec ( aa0, gg0, b0, q0, c0, r0, tol ) 
%HINFBISEC Bisection method for solving Hinf-polynomial equations.
%
%     [LAM,F,N,S] = HINFBISEC ( A, G, B, Q, C, R, TOL )
%
%   System of equations:
%                S*S~ = LAM^2*A - G      (1)
%       B*R*N + F*S*Q = Fs*C
%   where Fs=F~, R is a special factor and deg(F) < deg(B).
%
%   Inputs:
%     A, G - symmetric quasipolynomials
%     B, Q - stable polynomials
%     C    - polynomial
%     TOL  - tolerance (default tol = 1e-6)
%
%   Outputs:
%     LAM  - optimal value of LAM 
%     F, N, S - full-order solution of system (1)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
    if ~exist('tol','var'), tol = 1e-6; end;
    if ~exist('r0','var'),  r0 = 1; end;
    [aa0,gg0,b0,q0,c0,r0] = compat ( aa0, gg0, b0, q0, c0, r0 );
    polnMode = isa(aa0,'poln'); 
    aa = aa0; gg = gg0; q = q0; c = c0; r = r0;
    if polnMode
      var = aa0.var; 
      aa = aa.coef; gg = gg.coef;
      b = b0*r;
      b = b.coef; q = q.coef; c = c.coef;
    else
      b = conv(b0, r);
    end;           
%------------------------------------------------------
%   Degrees of Sigma, P and Q (full order)
%------------------------------------------------------
    degSi = max(deg(aa)/2, deg(gg)/2);
    degP  = deg(b0);
    degAll = degP + max(degSi+deg(q),deg(c));
    degQ  = degAll - deg(b);    
%------------------------------------------------------
%   Number of equations for full order solution
%------------------------------------------------------
%    n1 = delA + 1;
    n2 = degAll + 1;    
%------------------------------------------------------
%   Number of coefs for reduced order solution
%------------------------------------------------------
%    mSigma = delA + 1;  % deg(sigma) = delA  
    mP     = degP + 1;  % deg(P) = delB
    mQ     = degQ + 1;  % deg(Q) = delC
%    mQP = mQ + mP;    
%------------------------------------------------------
%   Find minimal lam as Hinf-norm of gam/a
%------------------------------------------------------
    diff = (deg(gg) - deg(aa)) / 2;
    if diff > 0
         R = tf ( gg, [aa zeros(1,diff)], 1 )';
    else R = tf ( [gg zeros(1,-diff)], aa, 1 )'; 
    end;
    xi0 = 1 / sqrt(zinfnorm(R));
%------------------------------------------------------
%   Check non-generic case
%------------------------------------------------------
    ag2 = sumpol2(aa,  -xi0^2*gg);
    if sum(ag2) < 1e-8*norm(ag2)
         sigma0 = sfactor ( ag2, 'd' );        
    else sigma0 = sfactfft ( ag2, 'd' ); 
    end;
    [P,Q] = dioph2 ( conv(sigma0,q), -xi0*c, b, degP );  
    if all(abs(roots(P)) < 1)
       [P1,Q1,err] = dioph2 ( -conv(sigma0,q), -xi0*c, b, degP );  
       lam = 1/xi0; count = 0;
       sigma = lam*sigma0;
       Q = lam*Q;
       if polnMode
         P = poln(P, var); P1 = poln(P1, var); 
         Q = poln(Q, var);  
         sigma = poln(sigma, var);  
       end;
       P = {P P1};
       return;
    end;
%------------------------------------------------------
%   Indices of unknown coeffieients in total vector
%------------------------------------------------------
    indQ = 1:mQ;
    indP = mQ+1:mQ+mP-1;        
%#######################################################
%   MAIN LOOOP: binary search
%#######################################################
    xiMin = 0;
    xiMax = xi0;
    count = 0;
    while xiMax - xiMin > tol
        count = count + 1; 
        xi = (xiMin + xiMax) / 2; 
        %lam = 1/xi
      %------------------------------------------------------
      %   Find sigma by spectral factorization
      %------------------------------------------------------
        sigma = sfactfft ( sumpol2(aa,  - xi*xi*gg), 'd' );        
      %------------------------------------------------------
      %   Construct A and B matrices for system of linear equations
      %------------------------------------------------------
        JQ = toep(b, n2, mQ);
        JP = toep(conv(sigma,q), n2, mP) - hank(xi*c, n2, mP);
        A  = [JQ JP];
        B  = - A(:,end);
        A  = A(:,1:end-1);
      %------------------------------------------------------
      %   Solve system of linear equations
      %------------------------------------------------------
        [flag,X] = linsys ( A, B, 'qr', 'refine' );
        if ~flag, X = A \ B; end;
      %------------------------------------------------------
      %   Extract polynomials Q and P
      %------------------------------------------------------
        Q = fliplr( X(indQ).' );
        P = fliplr( [X(indP); 1].' );    
      %-----------------JP-------------------------------------
      %   Bisection rule
      %------------------------------------------------------
        %disp([1/xi max(abs(roots(P)))])
        %disp(P)
        if any(abs(roots(P)) > 1), 
             xiMax = xi; 
        else xiMin = xi; 
        end;
    end        
%#######################################################
%   END OF MAIN LOOP
%#######################################################
    
%------------------------------------------------------
%   Extract full-order polynomials and xi
%------------------------------------------------------
    lam = 1/xi;
    Q = lam*Q;
    sigma = lam*sigma;
%------------------------------------------------------
%   Find reduce-order solution
%------------------------------------------------------
    rtsP = roots(P);
    rtsQ = roots(Q);
    err1 = min( abs(rtsP - 1) );
    err2 = min( abs(rtsP + 1) );
    rts0 = (err1 < err2) - (err1 > err2);
    errQ = min( abs(rtsQ - rts0) );
    rtsP = remove(rtsP, rts0);
    rtsQ = remove(rtsQ, rts0);
%------------------------------------------------------
%   Transform to polynomial mode
%------------------------------------------------------
    if polnMode
      P0 = poln(P, var); Q0 = poln(Q, var);  
      P = P(1)*poln(rtsP, ['r' var]);  
      Q = Q(1)*poln(rtsQ, ['r' var]);  
      sigma = poln(sigma, var);  
    else
      P0 = P; Q0 = Q;  
      P = P(1)*poly(rtsP);  
      Q = Q(1)*poly(rtsQ);  
    end;
    lam = abs(lam);
%------------------------------------------------------
%   Find solution error (before cancellation)
%------------------------------------------------------
    if polnMode
      p1 = sigma*sigma' - lam^2*aa0 + gg0;
      p2 = b0*r0*Q0 + P0*sigma*q0 - recip(P0)*c0;
    else        
      p1 = sumpol2(conv(sigma,fliplr(sigma)), -lam^2*aa0, gg0 );
      p2 = sumpol(conv(b,Q0), conv(P0,conv(sigma,q0)), ...
            -conv(fliplr(P0),c0));
    end;
    err = norm(p1) + norm(p2);
%------------------------------------------------------
%   Check common divisor of P and Q at the unit circle
%------------------------------------------------------
    errP = min(err1, err2);
    if max(errP,errQ) > 0.1, 
       lam = Inf; 
       return; 
    end;
%------------------------------------------------------
%   Check strict instability of P
%------------------------------------------------------
    if any(abs(roots(P)) > 1)
       lam = Inf; 
       return; 
    end;
       
%------- End of HINFBISEC.M --------- KYuP ----------
    
