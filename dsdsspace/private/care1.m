function [X,Poles,err,X1,X2] = care1 ( A, B, Q, R, S, E )
%CARE1 Solve continuous algebraic Riccati equation
%
%     [X,POLES,ERR,X1,X2] = CARE1 ( A, B, Q, R, S, E )
%
%   Inputs:
%     A, B, Q, R, S, E - constant matrices 
%
%   Outputs:
%     X     - solution of continuous algebraic Riccati equation
%                                       -1
%              A'XE + E'XA - (E'XB + S)R  (B'XE + S') + Q = 0
%
%	  POLES - poles of the optimal closed-loop system
%     ERR   - solution error	
%     X1,X2 - matrices such that X = X2 / X1
%

% References:
%   [1] W.F. Arnold, and A.J. Laub, ``Generalized Eigenproblem
%	    Algorithms and Software for Algebraic Riccati Equations,''
%	    Proc. IEEE, 72(1984), 1746--1754.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%   Check initial data
%------------------------------------------------------
    [n,n1] = size(A);
    if n~=n1, error('Matrix A must be square'); end;
    if any(size(Q)~=n)
       error('Matrices A and Q must have the same dimensions'); 
    end;    
    [nB,m] = size(B);
    if nB ~= n, error('Matrices A and B must have the same bnumber of rows'); end;
    if any(size(R)~=m)
       error('Matrix R must have compatible dimensions'); 
    end;    
    n2 = 2*n;   
    if ~exist('E','var'), 
       E = eye(n); 
       WithE = 0; 
    else
       WithE = 1; 
    end;    
    if ~exist('S','var'), S = zeros(n,m); end;    
%------------------------------------------------------
%   Extended pencil
%------------------------------------------------------
    a = [A  zeros(n) B
        -Q    -A'   -S
         S'    B'    R];
    b = [ E       zeros(n)
         zeros(n)    E'   ];
%------------------------------------------------------
%   Compression step when R is singular
%------------------------------------------------------
    [q,xx] = qr(a(:,n2+1:n2+m));
    a = q(:,n2+m:-1:m+1)'*a(:,1:n2);
    b = q(1:n2,n2+m:-1:m+1)'*b;   
%------------------------------------------------------
% 	Eigenstructure solution using QZ-algorithm
%------------------------------------------------------
    [aa,bb,q,z] = qz (a, b);   
%------------------------------------------------------
%   Find eigenvalues of the matrix pencil
%------------------------------------------------------
    D = diag(aa) ./ diag(bb);
%------------------------------------------------------
%	Order eigenvalues using complex Givens rotation
%   so that all eigs in the left half-plane are placed first
%------------------------------------------------------
    while 1
       OK = 1;
       for i=1:2*n-1
         if (real(D(i)) > 0)  &&  (real(D(i+1)) < 0)
            a11 = aa(i,i); a12 = aa(i,i+1); a22 = aa(i+1,i+1);
            b11 = bb(i,i); b12 = bb(i,i+1); b22 = bb(i+1,i+1);
   
            G1 = cgivens( a11*b12-a12*b11, a11*b22-a22*b11 );
            aa(i:i+1,:) = G1 * aa(i:i+1,:);
            bb(i:i+1,:) = G1 * bb(i:i+1,:);
            q(i:i+1,:)  = G1 * q(i:i+1,:);
   
            G2 = cgivens( a12*b22-a22*b12, a11*b22-a22*b11 );
            aa(:,i:i+1) = aa(:,i:i+1) * G2;
            bb(:,i:i+1) = bb(:,i:i+1) * G2;
            z(:,i:i+1)  = z(:,i:i+1) * G2;
       
            t = D(i); D(i) = D(i+1); D(i+1) = t;
            OK = 0;
          end
        end
        if OK, break; end;
     end
%------------------------------------------------------
%   For a specified E: get orthonormal basis
%------------------------------------------------------
    if WithE
       zz = E*z(1:n,1:n);
       a = [     zz
            z(n+1:n2,1:n)];
       [q,xx] = qr(a);
       z = q(:,1:n);
    end   
%------------------------------------------------------
%   Find solution and check symmetry
%------------------------------------------------------
    X1 = z(1:n,1:n);
    X2 = z(n+1:n2,1:n);
    X12 = X1'*X2;
    err = X12 - X12';
    if err > sqrt(eps)
        warning('Solution can be inaccurate')
    end
%------------------------------------------------------
%	Closed-loop poles
%------------------------------------------------------
	Poles = D(1:n);
%------------------------------------------------------
%	Form final solution
%------------------------------------------------------
   [L,U,P] = lu ( X1 );
   if rcond(U) > eps,
      X = real ( ((X2/U)/L)*P );
      X = (X + X') / 2;
   else
      X = inf*ones(n);
      err = -1;
      return;
   end
%------------------------------------------------------
%	Residual of the Riccati solution:
%------------------------------------------------------
    if rcond(R) < eps
       err = -1;                % --- singular R
    else
	   res = A'*X*E + E'*X*A - (E'*X*B + S)*(R\(B'*X*E + S')) + Q;
       err = norm(res) / max(1,norm(X));
    end;

%------- End of CARE1.M --------- KYuP ----------
