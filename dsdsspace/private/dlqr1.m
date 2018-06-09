function [K,S,Poles] = dlqr1 (A, B, Q, R, N)
%
% DLQR1  Linear-quadratic discrete regulator design
%
%   	[K,S,POLES] = DLQR1 ( A, B, Q, R, N )  
%
%   Inputs:
%	A, B, Q, R, N - constant matrices in standard LQ-problem
%			(cross-term N may be omitted, assumed zero)
%
%   Outputs:
%       K     - optimal gain matrix
%	    S     - solution to Riccati equation
%	    POLES - poles of the closed-loop system
%
%  NB: similar to DLQR from Control Toolbox
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 21-Sep-2006 $
%------------------------------------------------------
% 	Check input data
%------------------------------------------------------
	if ~exist('N','var')
	   N = zeros ( size(Q,1), size(R,2) );
	end;
%------------------------------------------------------
% 	Enforce symmetry and check positivity
%------------------------------------------------------
	Q = (Q + Q')/2;
	R = (R + R')/2;
%------------------------------------------------------
% 	Solve Riccati equation
%------------------------------------------------------
	[S,Poles] = dare1 ( A, B, Q, R, N );
    K = (B'*S*B + R) \ (B'*S*A + N');

%------- End of DLQR1.M --------- KYuP ----------
