function [M,E,P11,P22] = intaxa ( A, X, A2, T )
%INTAXA returns the symmetric integral 
%           M = int_0^T exp(A'*t)*X*exp(A*t) dt
%
%       M = INTAXA ( A, X, T )
%
%   Inputs: 
%      A - square n x n matrix
%      X - n x n matrix
%      T - integration limit
%
%   Outputs:
%      M = int_0^T exp(A*t)*X*exp(A'*t) dt 
%
%   Special syntax: 
%  -----------------
%       M = INTAXA ( A1, X, A2, T )
%
%   Inputs: 
%      A1, A2 - square matrices, n1 x n1 and n2 x n2
%      X - n1 x n2 matrix
%      T - integration limit
%
%   Outputs:
%           M = int_0^T exp(A1*t)*X*exp(A2'*t) dt
%

% [1] T. Chen, B.A. Francis, Optimal Sampled-Data Control Systems, 
%     Springer-Verlag, Berlin Heidelberg New York, 1995. Page 334.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 30-Oct-2006 $
%------------------------------------------------------
%       Check syntax
%------------------------------------------------------
        if nargin == 3, 
           T = A2; A2 = A; same = 1;     
        elseif nargin == 4 
           same = 0;
        else error('Function INTAXA requires 3 or 4 arguments');
        end;
        [nA,n1] = size(A);        
        if nA ~= n1, error('Matrix A must be square'); end;
        [nA2,n2] = size(A2);        
        if nA2 ~= n2, error('Matrix A2 must be square'); end;
        if any(size(X) ~= [nA nA2]), 
           error('Matrix X must have compatible dimensions'); 
        end;
%------------------------------------------------------
%       Perform integration using matrix exponential
%------------------------------------------------------
        AXA = [-A X;zeros(nA2,nA) A2.'];
        E = expm(AXA*T);
        P12 = E(1:nA,nA+1:end);
        if same
           P22 = E(nA+1:end,nA+1:end);
           M = P22.' * P12;
        else
           %P11 = E(1:nA,1:nA);
           %M = P11 \ P12;
           M = expm(A*T) * P12;
        end;
%------------------------------------------------------
%       Find P11 and P22 if required
%------------------------------------------------------
        if nargout > 2, P11 = E(1:nA,1:nA); end;
        if nargout > 3, P22 = E(nA+1:end,nA+1:end); end;

%------- End of INTAXA.M --------- KYuP ---------------
