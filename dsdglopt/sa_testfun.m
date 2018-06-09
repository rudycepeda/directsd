function y = sa_testfun ( x )
%SA_TESTFUN Complex function with many local minima.
%           Has a global minimum f(x,y)=0 at (0,0):
%
%    f(x,y)= a*x^2 + b*y^2 + c*[1-cos(alpha*x)] + d*[1-cos(gamma*y)] 
%
%   Parameters: a = 1; b = 2; 
%               c = 0.3; d = 0.4; 
%               alpha = 3*pi; gamma = 4*pi;
%
%   See also SIMANNEAL, SASIMPLEX.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        a = 1; b = 2; 
        c = 0.3; d = 0.4; 
        alpha = 3*pi; gamma = 4*pi;
        y = x(2); x = x(1);
        y = a*x^2 + b*y^2 + c*(1-cos(alpha*x)) + d*(1-cos(gamma*y));
    
%------- End of SA_TESTFUN.M --------- KYuP ----------

