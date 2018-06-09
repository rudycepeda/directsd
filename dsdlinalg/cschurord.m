function [u,t] = cschurord ( u, t, d )
%CSCHURORD Order eigenvalues in complex Schur form.
%          In increasing order of an auxiliary array D.
%
%     [U,T] = CSCHURORD ( U, T, D )
%
%   Inputs:
%     U - Hermitian transformaion matrix
%     T - complex Schur form
%     D - auxiliary array which determine the order
%
%   Outputs:
%     U - Hermitian transformaion matrix
%     T - complex Schur form
%
%   See also SCHUR, SCHURORD, CSF2RSF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        done = 0;
        n = size(t,1);
	    v = eye(n);
%------------------------------------------------------
%	Order eigenvalues in Schur decompisition
%   using complex Givens rotation
%------------------------------------------------------
        while ~done
           done = 1;     	   
           for k = 1:n-1
             if d(k) > d(k+1)
                done = 0;
            	g = cgivens ( t(k,k+1), t(k+1,k+1)-t(k,k) );
            	t(k:k+1,k:n) = g * t(k:k+1,k:n);
            	t(1:k+1,k:k+1) = t(1:k+1,k:k+1) * g';
            	v(:,k:k+1) = v(:,k:k+1) * g';
                t(k+1,k) = 0;
            	temp = d(k); d(k) = d(k+1); d(k+1) = temp;
             end
      	   end
   	    end
   	    u = u * v;

%------- End of CSCHURORD.M --------- KYuP ----------
