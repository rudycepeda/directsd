function [c,s] = cgivens (a, b)
%CGIVENS Complex Givens rotation matrix.
%        For some complex a and b find a real 'c' and
%        a complex 's' such that
%           G * | a | = | * |
%               | b |   | 0 |
%        and G is Hermitian (G~G = GG~ = I)
%               |  c    s |
%           G = | -s'   c |  
%
%     G     = CGIVENS ( A, B )
%     [C,S] = CGIVENS ( A, B )
%
%   Inputs:
%     A, B - rotation parameters
%
%   Outputs:
%     G   - rotation matrix
%     C,S - elements of rotation matrix
%
%   See also GIVENS.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Find c and s
%------------------------------------------------------
        if a == 0
           if b == 0
                c = 1; s = 0;
           else c = 0; s = 1; 
           end;
        else
           c = abs(a) / sqrt( abs(a'*a + b'*b) );
           s = c * b' / a';
        end
%------------------------------------------------------
%       Find G for special syntax
%------------------------------------------------------
        if nargout < 2
           G = [c s; -s' c];
           c = G;
        end;
        
%------- End of CGIVENS.M --------- KYuP ----------


