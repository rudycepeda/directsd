function H = hank ( P, R, C ) 
%HANK Hankel matrix for a polynomial.
%
%     H = HANK ( P, R, C )
%
%   Inputs:
%     P(s) - polynomial given by coefs [p_n p_(n-1) ... p1 p0]
%     R, C - row and column number of the matrix
%
%   Outputs:
%     H  - Hankel matrix:
%                  | 0 ...   0  p0 |
%              H = | 0 ...  p0  p1 |    
%                  | . ... ... ... |    
%                  | p0 .  ... pn  |    
%                  | 0 ...  pn  0  |    
%                  | 0 ...   0  0  |    
%                  | pn .    0  0  |    
%                  | 0 ...   0  0  |
%
%   See also TOEP, VEC.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        T = toep ( P, R, C );
        H = fliplr ( T );

%------- End of HANK.M --------- KYuP ----------
    
    
