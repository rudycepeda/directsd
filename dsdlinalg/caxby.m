function  [X,E,flag] = caxby ( A, B, Y )
%CAXBY Solution to liner matrix equation A*X*B = Y.
%       
%          [X,E,F] = CAXBY ( A, B, Y )
%
%  Inputs:
%       A, B, Y - constant matrices of compatible dimensions
%
%  Outputs:
%       X - constant solution to matrix equation A*X*B = Y or, 
%           if no solution exists, optimal approximant with minimal
%           Frobenius norm 
%       E - solution error
%       F - solvability flag (1 if an exact solution exists, 0 otherwise)
%
%   See also LINSYS.

% References:
%  [1] R.E.Skelton, T. Iwasaki, K.M. Grigoriadis, A unified approach
%      to linear control design, London, Taylor and Francis, 1998.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check solvability: equation has a solution if
%         1) A*A#*Y*B*B# = Y, or, equivalently,
%         2) (I - A*A#)*Y = 0, and Y*(I - B#*B) = 0
%       where A# and B# are Moore-Penrose pseudo-inverses
%
%       As a special case, linear equation A*x=y is solvable iff
%          (I - A*A#)*y = 0.
%       The set of all solutions is x = A#*y + (I - A*A#)*z, where
%       z is arbitrary vector
%------------------------------------------------------       
        Ai = pinv(A);
        Bi = pinv(B);
        flag = norm(A*Ai*Y*Bi*B - Y) < sqrt(eps);
%------------------------------------------------------
%       Find a solution
%          X = A#*Y*B# + Z - A#*A*Z*B*B#, where
%       Z is arbitrary matrix of compatible dimensions (here Z=0)
%
%       As a special case, the set of all solutions if linear equation A*x=y is 
%          x = A#*y + (I - A*A#)*z, where
%       z is arbitrary vector
%------------------------------------------------------       
        X = Ai * Y * Bi;
        E = A*X*B - Y;

%------- End of CAXBY.M --------- KYuP ---------------
