function M = inteat ( A, T )
%
% INTAXA returns the integral 
%           M = int_0^T exp(A*t) dt
%
%       M = INTEAT ( A, T )
%
%   Inputs: 
%      A - square n x n matrix
%      T - integration limit
%
%   Outputs:
%      M = int_0^T exp(A*t) dt 
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 30-Oct-2006 $
%------------------------------------------------------

	n = size(A,1);
	M = intaxa ( A, eye(n), zeros(n), T );

%------- End of INTEAT.M --------- KYuP ---------------
