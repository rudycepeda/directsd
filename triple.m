function [A,B,C,G] = triple ( A, B, C, tol )
%TRIPLE Extract greatest common divisor of there polynomials.
%
%     [A1,B1,C1,G] = TRIPLE ( A, B, C )
%     [A1,B1,C1,G] = TRIPLE ( A, B, C, TOL )
%                                             
%   Inputs:
%     A, B, C - coef arrays or polynomials
%     TOL  - tolerance (default 1e-8)
%
%   Outputs:
%     A1, B1, C1 - coprime triple
%     G - monic common polynomial factor
%
%   See also GCD, COPRIME.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Strip leading zeros
%------------------------------------------------------
        if ~exist('tol','var'), tol = sqrt(eps); end;         
        %[A,B,C] = compat ( A, B, C );
        polnMode = isa(A,'poln');
%------------------------------------------------------
%       Degenerate cases
%------------------------------------------------------              
        G = 1;
        if polnMode, G = poln(G,A.var); end;
        if deg(A) == 0  
           if norm(A) < tol, [B,C,G] = coprime ( B, C ); end;
           return;
        elseif deg(B) == 0
           if norm(B) < tol, [A,C,G] = coprime ( A, C ); end;
           return;
        elseif deg(C) == 0  
           if norm(C) < tol, [A,B,G] = coprime ( A, B ); end;
           return;
        end
%------------------------------------------------------
%       General case
%------------------------------------------------------              
        [A0,za] = delzero(A);
        [B0,zb] = delzero(B);
        [C0,zc] = delzero(C);
        if polnMode, 
           rtsA = A0.z; kA = A0.k;
           rtsB = B0.z; kB = B0.k;
           rtsC = C0.z; kC = C0.k;
        else
           rtsA = roots(A0); kA = A0(1);
           rtsB = roots(B0); kB = B0(1);
           rtsC = roots(C0); kC = C0(1);
        end;
        rtsG = [];
%------------------------------------------------------
%       Extract common roots
%------------------------------------------------------              
        i = 1;
        modified = 0;
        while i <= length(rtsA)
          R = rtsA(i); 
          tolR = max(tol, tol*abs(R));
          if length(rtsB) < 1  ||  length(rtsC) < 1, break; end;  
          errB = sort(abs(rtsB - R));  
          errC = sort(abs(rtsC - R));  
          if errB(1) < tolR  &&  errC(1) < tolR
            if abs(imag(R)) > eps % complex conjugate pair  
               R = [R; conj(R)];
            end;
            modified = 1;
            rtsG = [rtsG; R];
            rtsC = remove ( rtsC, R, tolR );
            rtsB = remove ( rtsB, R, tolR );
            rtsA = remove ( rtsA, R, tolR );
          else
            i = i + 1;
          end;
        end
%------------------------------------------------------
%       Final results
%------------------------------------------------------
        z0 = min([za zb zc]);
        za = za - z0;
        zb = zb - z0;
        zc = zc - z0;
        if modified || z0 > 0
          if polnMode
            u = poln([1 0], A.var);  
            A = kA*poln(rtsA, ['r' A.var])*u^za;
            B = kB*poln(rtsB, ['r' A.var])*u^zb;
            C = kC*poln(rtsC, ['r' A.var])*u^zc;
            G = poln(rtsG, ['r' A.var])*u^z0;
          else
            A = kA * [poly(rtsA) zeros(1,za)];
            B = kB * [poly(rtsB) zeros(1,zb)];
            C = kC * [poly(rtsC) zeros(1,zc)];
            G = [poly(rtsG) zeros(1,z0)];
          end;
        end;
        
%------- End of TRIPLE.M --------- KYuP ----------

