function [gamL,f] = dinfnorm(sys, tol)
%DINFNORM Compute Hinf-norm of the discrete-time model.
%
%     [N,F] = DINFNORM ( SYS, TOL )
%
%   Inputs:
%     SYS - discrete-time LTI model
%	  TOL - tolerance (default 1e-3)
%
%   Outputs:
%	  N - Hinf-norm of the closed-loop system
%	  F - frequency of the maximum gain
%
%   See also SDAHINORM, SDH2NORM.

% References:
%  [1] Bruisma, N.A., and M. Steinbuch, ``A Fast Algorithm to Compute
%      the Hinf-Norm of a Transfer Function Matrix,'' Syst. Contr. 
%      Letters, vol. 14, pp. 287-293, 1990.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1e-3; end;
        if ~isnumeric(tol), error('Incorrect TOL parameter'); end;
        if ~isdt(sys),
           error('SYS must be a discrete-time model'); 
        end;
%------------------------------------------------------
%       Special case: static gain
%------------------------------------------------------
        sys = minreal(ss(sys), sqrt(eps), 0);
        [a,b,c,d] = ssdata(sys);
        [xx,xx,K] = zpkdata(sys);
        if norm(K)<eps  ||  size(a,1) == 0,
           gamL = norm(d); 
           f = 0;
           return
        end
%------------------------------------------------------
%       Special case: infinite norm
%------------------------------------------------------
        eigs = eig(sys.a);
        [dist,i] = min(abs(1-abs(eigs)));
        if dist < sqrt(eps),
           gamL = Inf;  
           f = abs(angle(eigs(i)));
           return
        end
%------------------------------------------------------
%       Reduce to Hessenberg form
%------------------------------------------------------
        n = size(a,1);
        [u,a] = hess(a);
        b = u' * b;    
        c = c * u;
%------------------------------------------------------
%       First vector of test frequencies for modes
%------------------------------------------------------
        eigs_s = log(eigs(eigs~=0));
        w = abs(eigs_s);        
        ind = find(imag(eigs_s) >= 0 & w > eps);
        beta = abs(real(eigs)) ./ w;
        vf = w(ind) .* sqrt(max(0.25,1-2*beta(ind).^2));
        vfz = [exp(sqrt(-1)*vf) ; -1 ; 1];
%------------------------------------------------------
%       Compute lower bound
%------------------------------------------------------
        gamL = 0;
        for fz = vfz.',
          g = norm(d + c*((fz*eye(n) - a)\b));
          if g > gamL, 
             gamL = g;  
             f = abs(angle(fz)); 
          end
        end
%------------------------------------------------------
%       Gamma-iteration
%------------------------------------------------------
        relTol = tol;
        while 1
           g = (1+relTol) * gamL;         % upper bound candidate
           eigs = sympeig (a, b, c, d, g);
           mag = abs(eigs);
          %---------------------------------------------
          % Find eigenvalues at the unit circle and 
          % the corresponding frequencies
          %---------------------------------------------
           eig0 = eigs(abs(1-mag) < sqrt(eps));
           phi = angle(eig0);
           phi = phi(phi > 0);
          %---------------------------------------------
          % General step
          %---------------------------------------------
           if isempty(phi), 
              if abs(g - gamL) < tol, break; end;
              relTol = relTol / 10;
              continue; 
           end;
          %---------------------------------------------
          % New vector of test frequencies
          %---------------------------------------------
           gamL0 = gamL;
           len = length(phi);
           vfz = exp(sqrt(-1)*(phi(1:len-1)+phi(2:len))/2);
           for fz = vfz.',
             g = norm(d + c*((fz*eye(n) - a)\b));
             if g > gamL, 
                gamL = g;  
                f = abs(angle(fz)); 
             end
           end
          %---------------------------------------------
          % Check against infinite loop (undetected
          % eigs on the unit circle)
          %---------------------------------------------
           if gamL < gamL0*(1+sqrt(eps)), return; end;
        end
        
%------- End of DINFNORM.M --------- KYuP ----------

%######################################################
function eigs = sympeig (a, b, c, d, gamma)
%SYMPEIG Compute eigenvalues of the symplectic pencil.
%        Hinf-analysis of discrete-time systems. 
% References:
%   [1] T. Chen, B. Francis, Optimal Sampled-Data Control Systems, 
%       Springer-Verlag, Berlin etc., 1995.
%------------------------------------------------------
        n = size(a);
        c = c / gamma;    
        d = d / gamma;    
        [nout,nin] = size(d);
        R = eye(nout) - d*d';
        S = eye(nin) - d'*d;
        Ax = a + b*d'*(R\c);
        Sl = [   Ax    zeros(n)
              -c'*(R\c) eye(n)];
        Sr = [eye(n) -b*(S\b')
              zeros(n) Ax'];
        eigs = eig(Sl,Sr);

%------- End of SYMPEIG.M --------- KYuP ----------
