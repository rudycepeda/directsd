function [A11,E22,Q,Z] = weierstr ( A, E )
%WEIERSTR Canonical Weierstrass form for a matrix pencil sE - A.
%
%     [A1,E1,Q,Z] = WEIERSTR ( A, E )
%
%   Inputs: 
%     A, E - matrices of the same dimensions, E can be singular.
%
%   Outputs:
%     A1, E1 - matrices of the canonical Kroneker form
%     Q, Z   - transformation matrices such that
%                Q*(s*E-A)*Z = diag(s*I-A1, s*E1-I)  
%	           with E1 nilpotent so that E1^k=0 for some k>0		
%
%   See also KRON, QZPENCIL.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check input data
%------------------------------------------------------
        [n,n1] = size(A);
        if n ~= n1, error('Matrix A must be square'); end;
        if any(size(E) ~= n), error('Matrices A and E must have the same dimensions'); end;
%------------------------------------------------------
%       Separate infinite eigenvalues
%------------------------------------------------------
        [aa,ee,Q,Z] = qzpencil ( A, E, 'z' );       
        tol = 1.e-10;
        n1 = length( find(abs(diag(ee)) > tol) );
        n2 = n - n1;        
        A11 = aa(1:n1,1:n1);
        A12 = aa(1:n1, (n1+1):n);
        A22 = aa((n1+1):n, (n1+1):n);
        E11 = ee(1:n1,1:n1);
        E12 = ee(1:n1, (n1+1):n);
        E22 = ee((n1+1):n, (n1+1):n);
%------------------------------------------------------
%       Zeroing (1,2)-blocks using Sylvester equations
%         | I -L | | A11 A12 | | I R |   | A11 A11*R-L*A22+A12 | 
%         | 0  I | |  0  A22 | | 0 I | = |  0        A22       |
%       so that
%           A11*R - L*A22 = - A12
%           E11*R - L*E22 = - E12
%       by reduction to linear matrix equation
%  | kron(In2,A11) kron(-A22',In1) | |col(R)| = | -col(A12) |
%  | kron(In2,E11) kron(-E22',In1) | |col(L)|   | -col(E12) |
%------------------------------------------------------
        In1 = eye(n1);
        In2 = eye(n2);
        Z21 = zeros(n2,n1);

        M = [ kron(In2,A11) kron(-A22',In1)
              kron(In2,E11) kron(-E22',In1) ];
        n12 = numel(A12);                        
        b = - [ reshape(A12,n12,1)
                reshape(E12,n12,1)];
        x = M \ b;

        R = reshape( x(1:n12), n1, n2 );
        L = reshape( x((n12+1):end), n1, n2 );
        
        Q = [  In1 -L
               Z21 In2 ] * Q;
        Z = Z * [ In1  R
                  Z21 In2 ];        
%------------------------------------------------------
%       Transform A22 and E11 to identity matrices
%------------------------------------------------------
        Q1 = [ inv(E11) Z21'
                 Z21    In2 ];
        A11 = E11 \ A11; 
        Q = Q1 * Q; 
        E11 = eye(size(E11)); 
        
        Q2 = [ In1   Z21'
               Z21 inv(A22) ];
        E22 = A22\E22; 
        Q = Q2 * Q; 
        A22 = eye(size(A22));
%------------------------------------------------------
%       Transform matrices to real QZ-form
%------------------------------------------------------
        AA = [A11 Z21'
              Z21 A22];
        EE = [E11 Z21'
              Z21 E22];
        [QQ,Z,AA,EE] = csf2rsf(Q',Z,AA,EE);
        Q = QQ';
%------------------------------------------------------
%       Extract final result
%------------------------------------------------------   
        A11 = AA ( 1:n1, 1:n1 );
        E22 = EE ( (n1+1):n, (n1+1):n );
        for i=1:n2, E22(i,i) = 0; end;

%------- End of WEIERSTR.M --------- KYuP ----------           

