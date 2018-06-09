function [lam,P,Q,sigma,count,err] = hinfred ( aa0, gg0, b0, q0, c0, r0, tol ) 
%HINFRED Reduced-order solution for Hinf-polynomial equations.
%
%     [LAM,F,N,S] = HINFRED ( A, G, B, Q, C, R, TOL )
%
%   System of equations:
%              S*S~ = LAM^2*A - G      (1)
%       B*R*N + F*S*Q = Fs*C
%   where Fs=F~, R is a special factor and deg(F) < deg(B)
%
%   Inputs:
%     A, G - symmetric quasipolynomials
%     B, Q - stable polynomials
%     C, R - polynomials
%     TOL  - relative tolerance (default tol = 1e-4)
%
%   Outputs:
%     LAM  - optimal value of LAM 
%     F, N, S  - optimal solution of system (1)
%

% References: 
%  [1] Kwakernaak H., The polynomial approach to H-inf regulation // in
%      Lecture Notes in Mathematics, H-inf control theory, vol. 1496. 
%      London: Springer-Verlag, 1990. pp.~141-221.
%  [2] Saeki M., Method of solving a polynomial equation for an H-inf 
%      optimal control problem // IEEE Trans. Automat. Contr., vol. AC-34,
%      pp. 166-168, 1989.
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
    delA = max(deg(aa)/2, deg(gg)/2);
    delB = deg(b0);
    delC = max(delA+deg(q), deg(c));    
%------------------------------------------------------
%   Number of equations for reduced order solution
%------------------------------------------------------
    n1 = delA + 1;
    n2 = deg(b0) + delC;    
%------------------------------------------------------
%   Number of coefs for reduced order solution
%------------------------------------------------------
    mSigma = delA + 1;      % deg(sigma) = delA  
    mP     = delB;          % deg(P) = delB - 1
    mQ     = delC - deg(r); % deg(Q) = delC - 1
    mQP = mQ + mP;    
%------------------------------------------------------
%   Coefficient vectors with positive powers for
%           a*a~  and gamma*gamma~
%------------------------------------------------------
    Ca2 = aa(1:deg(aa)/2+1).';
    Ca2 = [zeros(n1-length(Ca2),1); Ca2];
    Cgam2 = gg(1:deg(gg)/2+1).';
    Cgam2 = [zeros(n1-length(Cgam2),1); Cgam2];    
%------------------------------------------------------
%   Find Sigma for xi=0 and form vec(Sigma)
%------------------------------------------------------
    sigma = sfactfft ( aa, 'd' );
    vecSigma = fliplr(sigma).'; 
    vecSigma = [vecSigma; zeros(mSigma-length(vecSigma),1)];    
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
    [P,Q] = dioph2 ( conv(sigma0,q), -xi0*c, b, mP );  
    if all(abs(roots(P)) < 1)
       [P1,Q1,err] = dioph2 ( -conv(sigma0,q), -xi0*c, b, mP );  
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
%   Find initial guess using generalized eigenvalue
%   approach by Saeki
%------------------------------------------------------
    TP  = toep(conv(sigma,q), n2, mP);
    TQ  = toep(b, n2, mQ);
    TwP = fliplr(toep(c, n2, mP)); 
       % form matrices for matrix pencil
    TA = [TQ TP];
    TB = [zeros(n2,mQ) TwP];
       % compress rows of TA and change TB accordingly
    [W,TA] = compress(TA,'row');
    TB = W*TB;
       % cut the tall matrices to square ones 
    TA = TA(1:mQP,:);
    TB = TB(1:mQP,:);
       % find maximal absolute value of 
       % finite pencil eigenvalues   
    eigXi = eig(TA, TB);
    eigXi = eigXi(~isinf(eigXi));
    xi = min(abs(eigXi)); 
%------------------------------------------------------
%   Solve linear system for coefs of Q and P
%   with unity highest coef of P 
%------------------------------------------------------
    indQ0 = 1:mQ;
    indP0 = mQ+1:mQ+mP-1;
    if xi > xi0, xi = 0.95*xi0; end;
    A = TA - xi*TB;       
    vecQP = - A(:,1:mQP-1) \ A(:,mQP);
%------------------------------------------------------
%   DeStabilize P
%------------------------------------------------------
    P = fliplr([vecQP(indP0); 1]);
    rtsP = roots(P);
    ind = find (abs(rtsP) > 1);
    if ~isempty(ind),
      rtsP(ind) = 1 ./ rtsP(ind);
      vecP = vec(poly(rtsP));
      vecQP(indP0) = vecP(1:end-1);
    end;
%------------------------------------------------------
%   Scaling
%------------------------------------------------------
    kScaleQ = 1.; %norm(vecQP(indQ0));
    b = b * kScaleQ;
    vecQP(indQ0) = vecQP(indQ0) / kScaleQ;
%------------------------------------------------------
%   Form the starting vector and define
%   intervals for coefs of all polynomials
%------------------------------------------------------
    X = [xi; vecSigma; vecQP];
    indXi    = 1;
    indSigma = 2:mSigma+1;
    indQ     = mSigma+2:mSigma+1+mQ;
    indP     = mSigma+2+mQ:mSigma+mQ+mP;  
%#######################################################
%   MAIN LOOOP: Newton-Raphson iterations
%#######################################################
    oldXi = Inf;
    count = 0; 
    reset = 0;
    while abs(oldXi-xi) > tol*norm(xi)
        count = count + 1; 
        if count > 200, break; end;
        %disp(1/xi)
      %------------------------------------------------------
      %   Extract reduce-order polynomials and xi
      %------------------------------------------------------
        vecSigma = X(indSigma);
        vecQ  = X(indQ);
        vecP  = [X(indP); 1];
        sigma = fliplr(vecSigma.');
        xi    = X(indXi);
        %Q = fliplr(vecQ.');
        P = fliplr(vecP.');    
      %------------------------------------------------------
      %   Find values of function vector
      %------------------------------------------------------
        F1 = fliplr(toep(sigma,n1,n1))*vecSigma - Ca2 + xi^2*Cgam2;
        F2 = toep(b, n2, mQ)*vecQ + toep(conv(sigma,q),n2,mP)*vecP ...
             - fliplr(toep(xi*c, n2, mP))*vecP;
        F = [F1; F2];
      %------------------------------------------------------
      %   Find Jacobian with respect to X, 
      %   except for the highest coefficient of P  
      %------------------------------------------------------
        Jsigma1 = toep(fliplr(sigma),n1,n1) + hank(sigma,n1,n1);
        Jxi1 = 2*xi*Cgam2;
        Jsigma2 = toep(conv(P,q),n2,mSigma);
        Jxi2 = - hank(c, n2, mP)*vecP;
        JQ = toep(b, n2, mQ);
        JP = toep(conv(sigma,q), n2, mP) - hank(xi*c, n2, mP);
        JP = JP(:,1:end-1);
    
        J = [Jxi1 Jsigma1 zeros(n1,mQ+mP-1)
             Jxi2 Jsigma2 JQ   JP      ];         
      %------------------------------------------------------
      %   Form new approximation
      %------------------------------------------------------
        [flag,dX] = linsys ( J, F, 'svd', 'refine' );
        if ~flag, dX = J \ F; end;
        X = X - dX;
        oldXi = xi;
        xi = X(indXi);
      %------------------------------------------------------
      %   Reset if P ceased to be strictly antistable
      %------------------------------------------------------
        if any(abs(roots(P)) > 1)  &&  ~reset
           xi = 0.01*xi0; 
           A = TA - xi*TB;       
           vecQP = - A(:,1:mQP-1) \ A(:,mQP);
           sigma = sfactfft ( sumpol2(aa, -xi^2*gg), 'd' );
           vecSigma = fliplr(sigma).'; 
           X = [xi; vecSigma; vecQP];
           reset = 1; 
        end;
    end        
%#######################################################
%   END OF MAIN LOOP
%#######################################################
    
%------------------------------------------------------
%   Extract reduce-order polynomials and xi
%------------------------------------------------------
    lam = 1/xi;
    sigma = fliplr(X(indSigma).');
    b = b / kScaleQ;
    Q = fliplr(X(indQ).') * kScaleQ;
    P = fliplr([X(indP); 1].');
    Q = lam * Q;
    sigma = lam * sigma;
%------------------------------------------------------
%   Make lam positive
%------------------------------------------------------
    if lam < 0, lam = - lam; end;
%------------------------------------------------------
%   Transform to polynomial mode
%------------------------------------------------------
    if polnMode
      P = poln(P, var);  
      Q = poln(Q, var);  
      sigma = poln(sigma, var);  
    end;
%------------------------------------------------------
%   Find solution error
%------------------------------------------------------
    if polnMode
      p10 = lam^2*aa0 - gg0;
      p1 = sigma*sigma' - p10;
      p20 = recip(P)*c0;
      p2 = b0*r0*Q + P*sigma*q0 - p20;
    else        
      p10 = sumpol2(lam^2*aa0, -gg0);  
      p1 = sumpol2(conv(sigma,fliplr(sigma)), -p10 );
      p20 = conv(fliplr(P),c0);
      p2 = sumpol(conv(b,Q), conv(P,conv(sigma,q0)), -p20);
    end;
    err = norm(p1)/norm(p10) + norm(p2)/norm(p20);
%------------------------------------------------------
%   Check full-order solution: 
%     P and Q must have common root at the unit circle
%------------------------------------------------------
    if polnMode
       [Px,Qx] = dioph2 ( sigma*q0, -c0, b, mP );  
       rtsP = Px.z;
       rtsQ = Qx.z;
    else        
       [Px,Qx] = dioph2 ( conv(sigma,q0), -c0, b, mP );  
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
    if any(abs(roots(P)) > 1)
       lam = Inf; 
       return; 
    end;

%------- End of HINFRED.M --------- KYuP ----------
    
