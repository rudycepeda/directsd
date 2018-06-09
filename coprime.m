function [A,B,G] = coprime ( A, B, tol )
%COPRIME Extract greatest common divisor of two polynomials.
%
%     [AX,BX,G] = COPRIME ( A, B )
%     [AX,BX,G] = COPRIME ( A, B, TOL )
%
%   Inputs:
%     A, B - polynomials
%     TOL  - tolerance (default 1e-8)
%
%   Outputs:
%     AX, BX - coprime polyniomials
%     G      - monic common polynomial factor
%
%   See also GCD, TRIPLE.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%         Check data  
%------------------------------------------------------
        if ~exist('tol','var'), tol = sqrt(eps); end;         
        [A,B] = compat ( A, B );
        polnMode = isa(A,'poln');
        G = 1;
        if polnMode, G = poln(G,A.var); end;
        if deg(A) < 1  || norm(A) < 1e-10, return; end;
        if deg(B) < 1  || norm(B) < 1e-10, return; end;       
%------------------------------------------------------
%         Find roots
%------------------------------------------------------
        if polnMode, 
           rtsA = A.z; kA = A.k;
           rtsB = B.z; kB = B.k;
        else
           [A0,za] = delzero(A);
           [B0,zb] = delzero(B);
           rtsA = roots(A0); kA = A0(1);
           rtsB = roots(B0); kB = B0(1);
        end;
        rtsG = [];
%------------------------------------------------------
%         Find common roots
%------------------------------------------------------
        i = 1;
        modified = 0;
        while i <= length(rtsA)
          R = rtsA(i);  
          tolR = max(tol, tol*abs(R));
          if length(rtsB) < 1, break; end;  
          errB = sort(abs(rtsB - R));  
          if errB(1) < tolR
            if abs(imag(R)) > eps % complex conjugate pair  
               R = [R; conj(R)];
            end;
            rtsG = [rtsG; R]; 
            if tolR < 1e-10, keyboard; end;
            rtsA = remove ( rtsA, R, tolR );
            rtsB = remove ( rtsB, R, tolR );
            modified = 1;
          else   
            i = i + 1;
          end;
        end
%------------------------------------------------------
%       Final results
%------------------------------------------------------
        if modified
          if polnMode
            A = kA*poln(rtsA, ['r' A.var], A.shift);
            B = kB*poln(rtsB, ['r' A.var], B.shift);
            G = poln(rtsG, ['r' A.var]);
          else
            A = kA * [poly(rtsA) zeros(1,za-zb)];
            B = kB * [poly(rtsB) zeros(1,zb-za)];
            G = [poly(rtsG) zeros(1,min(zb,za))];
          end;
        end;

%------- End of COPRIME.M --------- KYuP ---------------
