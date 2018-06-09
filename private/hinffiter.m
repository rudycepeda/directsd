function [lam,F,N,sigma,count1,count2,err] = hinffiter ( aa0, gg0, b0, q0, c0, r0, tol ) 
%HINFFITER F-iteration method for solving Hinf-polynomial equations.
%
%     [LAM,F,N,S] = HINFFITER ( A, G, B, Q, C, R, TOL )
%
%   System of equations:
%                S*S~ = LAM^2*A - G      (1)
%       B*R*N + F*S*Q = Fs*C
%   where Fs=F~, R is a special factor and deg(F) < deg(B).
%
%   Inputs:
%     A, G - symmetric quasipolynomials
%     B, Q - stable polynomials
%     C, R - polynomials
%     TOL  - tolerance (default tol = 1e-6)
%
%   Outputs:
%     LAM  - optimal value of LAM 
%     F, N, S  - optimal solution of system (1)
%

% References: 
% [1] M. Grimble, Robust Industrial Control, UK: Prentice Hall,
%     Hemel Hempstead, 1994, p. 157.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        %aa0 = aa0.coef; g0 = g0.coef; b0 = b0.coef; 
        %q0 = q0.coef; c0 = c0.coef; r0 = r0.coef; 
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
        delA = max(deg(aa)/2, deg(gg)/2);
        delB = deg(b0);
        delC = max(delA+deg(q), deg(c));    
%------------------------------------------------------
%   Number of equations for reduced-order solution
%------------------------------------------------------
%    n1 = delA + 1;
        n2 = delB + delC;
        mF = delB;               % deg(F) = delB - 1
        mN = delC - deg(r);      % deg(N) = delC - 1 - degR
        indN = 1:mN;
        indF = mN+1:mN+mF;        
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
        [F,N] = dioph2 ( conv(sigma0,q), -xi0*c, b, mF );  
        if all(abs(roots(F)) < 1)
           [F1,N1,err] = dioph2 ( -conv(sigma0,q), -xi0*c, b, mF );  
           lam = 1/xi0; count1 = 0; count2 = 0;
           sigma = lam*sigma0;
           N = lam*N;
           if polnMode
             F = poln(F, var); F1 = poln(F1, var); 
             N = poln(N, var);  
             sigma = poln(sigma, var);  
           end;
           F = {F F1};
           return;
        end;
%------------------------------------------------------
%       Start from lam=Inf
%------------------------------------------------------
        lam = 1e8;
        Right = conv(c, tol);        
%------------------------------------------------------
%       Main iteration loop
%------------------------------------------------------
        count1 = 0; count2 = 0;
        while 1
          sigma = sfactfft ( sumpol2(aa,  - gg/lam^2), 'd' );
          %disp(lam), %pause
          lamOld = lam;
          while 1
             count2 = count2 + 1;
             if count2 > 20, break; end;
             lamOldX = lam; 
            %------------------------------------------------------
            %   Construct A and B matrices for system of linear equations
            %------------------------------------------------------
            %[Flam,N] = dioph ( conv(sigma,q), b, Right ); 
             JN = toep(b, n2, mN);
             JF = toep(conv(sigma,q), n2, mF);
             B  = toep(Right, n2, 1);
            %------------------------------------------------------
            %   Solve system of linear equations
            %------------------------------------------------------
             X  = [JN JF] \ B;
            %------------------------------------------------------
            %   Find F(z) and Fs(z)
            %------------------------------------------------------
             N = fliplr(X(indN).'); 
             Flam = fliplr(X(indF).'); 
             lam = abs(Flam(1));
             F = Flam / lam;
             Fs = fliplr( F );
             Right = conv(Fs, c);        
             if abs(lamOldX-lam) < tol*lam, break; end;
          end;
          if abs(lamOld-lam) < tol*lam, break; end;
          count1 = count1 + 1;
          if count1 > 200, break; end;
        end;  
%------------------------------------------------------
%       Transform back to lam
%------------------------------------------------------
        sigma = lam*sigma;
        if polnMode
          F = poln(F, var);  
          N = poln(N, var);  
          sigma = poln(sigma, var);  
        end;
        lam = abs(lam);
%------------------------------------------------------
%       Check error
%------------------------------------------------------
        if polnMode
           p10 = lam^2*aa0 - gg0;
           p1 = sigma*sigma' - p10;
           p20 = recip(F)*c0;
           p2 = b0*r0*N + F*sigma*q0 - p20;
        else        
           p10 = sumpol2(lam^2*aa0, -gg0);  
           p1 = sumpol2(conv(sigma,fliplr(sigma)), -p10 );  
           p20 = conv(fliplr(F),c0);
           p2 = sumpol(conv(b,N), conv(F,conv(sigma,q0)), -p20);
        end;
        err = norm(p1)/norm(p10) + norm(p2)/norm(p20);
%------------------------------------------------------
%   Check full-order solution: 
%     P and Q must have common root at the unit circle
%------------------------------------------------------
        if polnMode
           [Px,Qx] = dioph2 ( sigma*q0, -c0, b, mF );  
           rtsP = Px.z;
           rtsQ = Qx.z;
        else        
           [Px,Qx] = dioph2 ( conv(sigma,q0), -c0, b, mF );  
           rtsP = roots(Px);
           rtsQ = roots(Qx);
        end;
        err1 = min( abs(rtsP - 1) );
        err2 = min( abs(rtsP + 1) );
        errP = min(err1, err2);
        rts0 = (err1 < err2) - (err1 > err2);
        errQ = min( abs(rtsQ - rts0) );
        if max(errP,errQ) > 0.1, 
           lam = Inf; 
           return; 
        end;
%------------------------------------------------------
%   Check strict instability of P
%------------------------------------------------------
        if any(abs(roots(F)) > 1)
           lam = Inf; 
           return; 
        end;
    
%------- End of HINFFITER.M --------- KYuP ----------
    
    
