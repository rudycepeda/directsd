function [P,Poles,err,p1,p2] = care2 ( A, Q, R )
%CARE2 Solve continuous-time algebraic Riccati equation.
%
%     [P,POLES,ERR,P1,P2] = CARE2 ( A, Q, R )
%
%   Inputs:
%     A, Q, R - constant matrices 
%
%   Outputs:                                              
%     P     - solution of Riccati equation A'P + PA - PRP + Q = 0
%     POLES - poles of the optimal closed-loop system
%	  ERR   - solution error
%     P1,P2 - matrices such that P = P2 / P1
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
% 	Special case: A is stable, Q = 0 => P = 0
%------------------------------------------------------
	n = size(A);
	if isempty( find(eig(A)>0,1) ) &&  norm(Q) == 0
       p1 = eye(n);
       p2 = zeros(n);
       P  = zeros(n); 
       Poles = eig(A); 
       err = 0;
	   return
	end
%------------------------------------------------------
% 	Hamiltonian matrix
%------------------------------------------------------
	H = [A -R;-Q -A'];
%------------------------------------------------------
% 	Eigenstructure solution
%------------------------------------------------------
   	[v,dd] = reig(H);  	
   	p1 = v ( 1:n, 1:n );
   	p2 = v ( n+1:2*n, 1:n);
  	P = real( p2 / p1 );
%------------------------------------------------------
%	Closed-loop poles
%------------------------------------------------------
	Poles = eig(dd);
    [xx,ind] = sort( real(Poles) );
    Poles = Poles(ind(1:n));
%------------------------------------------------------
%	Residual of the Riccati solution:
%------------------------------------------------------
	err = P*A + A'*P - P*R*P + Q;

%------- End of CARE2.M --------- KYuP ----------
