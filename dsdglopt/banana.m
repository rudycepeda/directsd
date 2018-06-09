function y = banana ( x )
%BANANA Rosenbrock's "banana function" f(x,y)= 100*(y-x^2)^2+(1-x)^2.
%
%   This function has a unique minimum at the
%   point (x=1,y=1) where f(x,y)=0.
%
%   See also OPTGLOB.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------

	y = (100*(x(2) - x(1)*x(1))^2 + (1-x(1))^2);

%------- End of BANANA.M --------- KYuP ----------

