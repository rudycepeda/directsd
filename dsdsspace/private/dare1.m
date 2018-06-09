function [X,Poles,err,X1,X2] = dare1 ( A, B, Q, R, S, E )
%DARE1 Solve DARE with possible unit-circle eigenvalues.
%
%     [X,POLES,ERR,X1,X2] = DARE1 ( A, B, Q, R, S, E )
%
%   Inputs:
%     A, B, Q, R, S, E - constant matrices 
%
%   Outputs:
%     X - solution of discrete algebraic Riccati equation
%                                                    -1
%            A'*X*A - E'*X*E - (A'*X*B+S)(R + B'*X*B) *(B'*X*A+S') + Q = 0
%
%	  POLES - poles of the optimal closed-loop system
%     ERR   - solution error	
%     X1,X2 - matrices such that X = X2 / X1
%

% [1] W.F. Arnold, and A.J. Laub, ``Generalized Eigenproblem
%	  Algorithms and Software for Algebraic Riccati Equations,''
%	  Proc. IEEE, 72(1984), 1746--1754.
% [2] K. Polyakov, Separation of eigenvalues on the unit circle, 2005.  
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
    if nB ~= n, error('Matrices A and B must have the same number of rows'); end;
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
%   Scaling weight functions
%------------------------------------------------------   
    scale = norm(Q,1)+norm(R,1)+norm(S,1);
    QQ = Q / scale;   
    RR = R / scale;   
    SS = S / scale;
%------------------------------------------------------
%   Extended pencil
%------------------------------------------------------
    a = [  E      zeros(n) zeros(n,m)
         zeros(n)    A'    zeros(n,m)
         zeros(m,n) -B'    zeros(m)];
    b = [ A  zeros(n)     B
         -QQ    E'      -SS
         SS' zeros(m,n)  RR ];
%------------------------------------------------------
%   Compression step when R is singular
%------------------------------------------------------
    [q,xx] = qr(b(:,n2+1:n2+m));
    a = q(:,n2+m:-1:m+1)' * a(:,1:n2);
    b = q(:,n2+m:-1:m+1)' * b(:,1:n2);
%------------------------------------------------------
% 	Eigenstructure solution using QZ-algorithm
%------------------------------------------------------
    [aa,bb,q,z] = qz (a, b);
%------------------------------------------------------
%   Find eigenvalues of the matrix pencil
%------------------------------------------------------
    alpha = diag(aa);
    beta = diag(bb);   
    D = [];
    for i=1:2*n
       if abs(beta(i)) < 1.e-10
            D = [D; Inf];
       else D = [D; alpha(i)/beta(i)]; 
       end;
    end;
%------------------------------------------------------
%	Split eigenvalues at the unit circle
%------------------------------------------------------
    ucPoles = D*0;
    dir = -1;
    for i=1:2*n
       if abs(D(i)-1) < sqrt(eps), 
          D(i) = 1 + dir*eps;
          aa(i,i) = D(i) * bb(i,i);
          dir = - dir;
          ucPoles(i) = 1;
       end;
    end;
%------------------------------------------------------
%	Order eigenvalues using complex Givens rotation
%   so that all eigs outside the unit disk are placed first
%------------------------------------------------------
    while 1
       OK = 1;
       for i=1:2*n-1
         if abs(D(i)) < abs(D(i+1))
            [aa,bb,q,z,D,ucPoles] = movedown ( aa, bb, q, z, D, ucPoles, i );  
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
%   Find such a distribution of the unit-circle poles
%   that X1 is nonsingular
%------------------------------------------------------
    X1 = z(1:n,1:n); condX1 = cond(X1);
    X2 = z(n+1:n2,1:n);
    if condX1 > 1e12
       ucno = find(ucPoles > 0);
       if ~WithE  &&  mod(length(ucno),2) == 0 
          [X1,X2] = permUcPoles ( aa, bb, q, z, D, ucPoles, ...
                           min(ucno), max(ucno), X1, X2 );   
       end;
    end;
%------------------------------------------------------
%   Find solution and check symmetry
%------------------------------------------------------
    X12 = X1'*X2;
    err = X12 - X12';
    if err > sqrt(eps)
        warning('Solution can be inaccurate')
    end
%------------------------------------------------------
%	Closed-loop poles
%------------------------------------------------------
	Poles = D((n+1):end);
%------------------------------------------------------
%	Form final solution
%------------------------------------------------------
    [L,U,P] = lu ( X1 );
    if rcond(U) > eps,
      X = real ( ((X2/U)/L)*P );
      X = scale * (X + X') / 2;
    else
      X = inf*ones(n);
      err = -1;
      return;
    end
%------------------------------------------------------
%	Residual of the Riccati solution:
%
%------------------------------------------------------
%    RBX = R + B'*X*B
%    eig(RBX)
%------------------------------------------------------
    RBX = R + B'*X*B;
    if cond(RBX) < 1/eps 
      res = A'*X*A - E'*X*E - (A'*X*B+S)/ RBX *(B'*X*A+S') + Q;
      err = norm(res) / max(1,norm(X));
    else
      warning('Matrix R+B''*X*B is singular.');   
      err = 0;
    end;    

%------- End of DARE1.M --------- KYuP ----------
%################################################
function [X1,X2] = permUcPoles ( aa, bb, q, z, D, ucPoles, nFrom, nTo, X1, X2 )
%------------------------------------------------------
%
% PERMUCPOLES permutate poles at the unit circle to
%             improve conditioning of X1
%------------------------------------------------------
%       K. Yu. Polyakov         10 Feb 2005    
%                               
%------------------------------------------------------
%       Check errors
%------------------------------------------------------
        n = size(X1,1);
        condX1 = cond(X1);
        p = (nTo - nFrom + 1) / 2;
        if p ~= floor(p), return; end;        
%------------------------------------------------------
%       Final step: two elements
%------------------------------------------------------
        if p == 1
           for i=1:2 
             X1a = z(1:n,1:n); condX1a = cond(X1a);
             if condX1a < condX1
               X1 = X1a; condX1 = condX1a;
               X2 = z(n+1:2*n,1:n);
             end;
             [aa,bb,q,z,D,ucPoles] = movedown ( aa, bb, q, z, D, ucPoles, nFrom );
           end;
           return;
        end;
%------------------------------------------------------
%       Loop for the left side
%------------------------------------------------------
        for i=1:p
          for j=1:p-1 
            [aa,bb,q,z,D,ucPoles] = movedown ( aa, bb, q, z, D, ucPoles, nFrom+j-1 );            
          end;
         %------------------------------------------------------
         %  Loop for the right side
         %------------------------------------------------------
          for ii=1:p
             for jj=1:p-1 
               [aa,bb,q,z,D,ucPoles] = moveup ( aa, bb, q, z, D, ucPoles, nTo-jj+1 );            
             end;
             [X1a,X2a] = permUcPoles ( aa, bb, q, z, D, ucPoles, nFrom+1, nTo-1, X1, X2 );
             condX1a = cond(X1a);
             if condX1a < condX1
               X1 = X1a; condX1 = condX1a; X2 = X2a;
             end;
          end;
        end;

%------- End of PERMUCPOLES.M --------- KYuP ----------

%################################################
function [aa,bb,q,z,D,ucPoles] = moveup ( aa, bb, q, z, D, ucPoles, i )
%------------------------------------------------------
%
% MOVEUP move i-th eigenvalue one position up
%          using complex Givens rotation
%------------------------------------------------------
%       K. Yu. Polyakov         10 Feb 2005    
%                               
%------------------------------------------------------
    [aa,bb,q,z,D,ucPoles] = movedown ( aa, bb, q, z, D, ucPoles, i-1 );
%------- End of MOVEUP.M --------- KYuP ----------
%################################################
function [aa,bb,q,z,D,ucPoles] = movedown ( aa, bb, q, z, D, ucPoles, i )
%------------------------------------------------------
%
% MOVEDOWN move i-th eigenvalue one position down
%          using complex Givens rotation
%------------------------------------------------------
%       K. Yu. Polyakov         10 Feb 2005    
%                               
%------------------------------------------------------
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
    t = ucPoles(i); ucPoles(i) = ucPoles(i+1); ucPoles(i+1) = t;

%------- End of MOVEDOWN.M --------- KYuP ----------
