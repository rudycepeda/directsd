function [X,Phi1,Phi2] = intaba ( A1, B, A2, t )
%INTABA Special integral from matrix exponentials.
%       These appear in state-space representation of
%       discrete models of continuous plants
%
%     [X,PHI1,PHI2] = INTABA (A1, B, A2, t)
%
%   Inputs:
%     A1, A2 - square constant matrices
%     B      - constant matrix suitable for A1*B*A2 product
%     t      - integration limit
%
%   Outputs:
%     X = exp(A1*t) * int[0..t] exp(-A1*v) * B * exp(A2*v) dv
%       = int[0..t] exp(A1*v) * B * exp(-A2*v) dv * exp(A2*t)
%     PHI1 = exp(A1*t)
%     PHI2 = exp(A2*t)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check correctness
%------------------------------------------------------
        [r1,c1] = size(A1);
        [r2,c2] = size(A2);
        [rB,cB] = size(B);
        if r1 ~= c1, error('Matrix A1 must be square'); end;
        if r2 ~= c2, error('Matrix A2 must be square'); end;
        if (c1 ~= rB)  ||  (cB ~= r2), error('Matrix B has incorrect dimensions'); end;   
%------------------------------------------------------
%       Form basic matrix
%------------------------------------------------------
        Z = [A1 B; zeros(r2,c1) A2];
        EZ = expm(Z*t);
        X = EZ(1:r1,r1+1:end);        
          % Phi1 = EZ(1:r1,1:r1);
          % Phi2 = EZ(r1+1:end,r1+1:end);
        Phi1 = expm(A1*t);
        Phi2 = expm(A2*t);

%------- End of INTABA.M --------- KYuP ----------           
