function [c,s,G] = givens (a, b)
%GIVENS Real Givens rotation matrix.
%       For some real a and b find real c and s
%       such that  G * | a | = | * |
%                      | b |   | 0 |
%       for the Givens rotation matrix 
%             |  c   s  |
%         G = | -s   c  |
%
%     G     = GIVENS ( A, B )
%     [C,S] = GIVENS ( A, B )
%
%   Inputs:
%     A, B - rotation parameters
%
%   Outputs:
%     G   - rotation matrix
%     C,S - elements of rotation matrix
%
%   See also CGIVENS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        if b == 0
           c = 1; s = 0;
        else
           if abs(b) > abs(a)
              tau = - a / b;
              s = 1 / sqrt(1 + tau*tau);
              c = s * tau;
           else
              tau = - b / a;
              c = 1 / sqrt(1 + tau*tau);
              s = c * tau;
           end;
        end;
%------------------------------------------------------
%       Return Givens matrix rather than c & s
%------------------------------------------------------
        G = [c -s; s c];
        if nargout <= 1
           c = G;     
        end;

%------- End of GIVENS.M --------- KYuP ----------           
