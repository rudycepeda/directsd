function [K,S,Poles] = lqr1 (A, B, Q, R, N)
%LQR1 Linear-quadratic regulator design.
%
%     [K,S,E] = LQR1 ( A, B, Q, R, N )  
%
%   Inputs:
%	  A, B, Q, R, N - constant matrices in standard LQ-problem
%			         (cross-term N may be omitted, assumed zero)
%   Outputs:
%     K     - optimal gain matrix
%	  S     - solution to Riccati equation
%	  POLES - poles of the closed-loop system
%
%  NB: Similar to LQR from Control Toolbox, but based on CARE1.
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
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
	[S,Poles] = care1 ( A, B, Q, R, N );
    K = R \ (N' + B'*S);

%------- End of LQR1.M --------- KYuP ----------

