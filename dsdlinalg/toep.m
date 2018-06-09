function T = toep ( a, r, c )
%TOEP  Toeplitz matrix for a polynomial.
%
%     T = TOEP ( A, R, C )
%
%   Inputs:
%     A    - a polynomial
%     R, C - row and column number of the desired matrix
%
%   Outputs:
%     T - Toeplitz matrix
%           | A(0)  0    0   ... |
%           | A(1) A(0)  0   ... |
%       T = |  ... A(1) A(0) ... |
%           |      ...  A(1) ... |
%           |           ...      |
%   If A contains two rows, say A1 and A2, then:
%           | A1(0)   0     0   ... |
%           | A1(1) A2(0)   0   ... |
%       T = |  ...  A2(1) A1(0) ... |
%           |       ...   A1(1) ... |
%           |              ...      |
%
%   See also HANK, VEC.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        ra = size ( a, 1 );
        if ra == 1,
           a1 = a;
           a2 = a;
        else
           a1 = a(1,:);
           a2 = a(2,:);
        end
        if isa(a,'poln'),
          a1 = a1.coef; a2 = a2.coef;  
        end;
        if nargin < 4, 
          a1 = fliplr(a1); 
          a2 = fliplr(a2); 
        end;
%------------------------------------------------------
%       Define sizes
%------------------------------------------------------
        n1 = length(a1);
        n2 = length(a2);
%------------------------------------------------------
%       First block
%------------------------------------------------------
        T = zeros ( r, c );
        if isa(a,'sym'), T = sym(T); end;
        for i = 1:2:c
           last = min(n1, r-i+1);
           if last > 0
              T(i:i+last-1,i) = a1(1:last).';
           end;
        end
%------------------------------------------------------
%       Second block
%------------------------------------------------------
        for i = 2:2:c
           last = min(n2, r-i+1);
           if last > 0
              T(i:i+last-1,i) = a2(1:last).';
           end;
        end        

%------- End of TOEP.M --------- KYuP ----------
