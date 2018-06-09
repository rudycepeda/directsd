function [A1,B1,Q,Z,D] = qzpencil ( A, B, type )
%QZPENCIL QZ-decomposition for a matrix pencil.
%
% Ordered real Schur transformation for a square matrix A
% --------------------------------------------------------
%
%     [A1,Q,D] = QZPENCIL ( A );      
%     [A1,Q,D] = QZPENCIL ( A, TYPE );
%
%   Inputs:  
%     A    - square matrix
%     TYPE - type of sorting of
%         's' - ascending real parts (default)
%         'i' - descending real parts
%         'z' - ascending absolute values
%         'd' - descending absolute values
%
%   Outputs: 
%     A1 - transformed upper block triangular matrix
%     Q  - unitary transformation matrix such that Q'*A*Q = A1
%     D  - reordered eigenvalues  
%
% QZ-transformation for a matrix pencil (A - lam*B)
% --------------------------------------------------------
%
%     [A1,B1,Q,Z,D] = QZPENCIL ( A, B );    %  no sorting
%     [A1,B1,Q,Z,D] = QZPENCIL ( A, B, TYPE );
%       
%   Inputs:  
%     A, B - square matrices of the same size
%     TYPE - type of sorting of
%         's' - ascending real parts (default)
%         'i' - descending real parts
%         'z' - ascending absolute values
%         'd' - descending absolute values
%       For type 'd', infinite eigenvalues are placed first,
%       otherwise infinite eigenvalues are placed last,
%	      's0' - the same as 's', but Inf eigs are in the center
%		  'i0' - the same as 'i', but Inf eigs are in the center
%		  'z0' - the same as 'z', but Inf eigs are in the center
%		  'd0' - the same as '0', but Inf eigs are in the center
%
%   Outputs: 
%     A1, B1 - transformed upper triangular matrices
%     Q, Z   - unitary transformation matrices such that
%                  Q*A*Z = A1 and Q*B*Z = B1  
%     D      - reordered eigenvalues
%
%   See also QZ, CSF2RSF.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check syntax 
%------------------------------------------------------
        [n,c] = size ( A );
        if n ~= c, error( 'Matrix A must be square' ); end;
        onlyA = 0;
        if nargin == 1
           type = 's'; B = eye(n); onlyA = 1;
        end;
        if ischar(B)
           if nargin > 2
              error('Incorrect number of input arguments');  
           end;     
           type = B; B = eye(n); onlyA = 1;
        end;       
        if onlyA  &&  nargout > 3
           error('Incorrect number of input arguments');
        end;
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        [r,b] = size ( B );
        if r ~= b, error( 'Matrix B must be square' ); end;
        if n ~= r
           error( 'Matrices A and B must have the same size' );
        end;
        if exist('type','var')
          if ~isequal(type,'s') && ~isequal(type,'i') && ...
             ~isequal(type,'z') && ~isequal(type,'d') && ...
             ~isequal(type,'s0') && ~isequal(type,'i0') && ...
             ~isequal(type,'z0') && ~isequal(type,'d0') 
             error('Incorrect sorting parameter');
          end;
        else type = 's'; 
        end;
        if length(type) < 2, 
          if type == 'd',
               type = [type '1'];       % infinite eigs first
          else type = [type '2'];       % infinite eigs last
          end;    
        end;
%------------------------------------------------------
%       Schur transformation for a single matrix 
%------------------------------------------------------
        if onlyA
           [Q,A1] = schur ( A );
           [Z,A1] = rsf2csf ( Q, A1 );
           B1 = B;
           Q = Z';
           D = diag(A1);
        else
%------------------------------------------------------      
%       QZ-transformation for matrix pencil
%------------------------------------------------------
           [A1,B1,Q,Z] = qz ( A, B );          
           qzcheck(A,A1,Q,Z);
           qzcheck(B,B1,Q,Z);
          %--------------------------------------------           
          % Find eigenvalues of the matrix pencil
          %--------------------------------------------           
           alpha = diag(A1);
           beta = diag(B1);
           D = [];
           for i=1:n
              if abs(beta(i)) < 1.e-10
                   D = [D; Inf];
              else D = [D; alpha(i)/beta(i)]; 
              end;
           end;          
        end;       
%------------------------------------------------------      
%       Form ordering array
%------------------------------------------------------      
        switch type(1)
          case 's', ordN = real(D);
          case 'i', ordN = - real(D);
          case 'z', ordN = abs(D) - 1;
          case 'd', ordN = - (abs(D) - 1);
        end;	
        ind = find ( isinf(D) );
        switch type(2)
          case '1', ordN(ind) = min(ordN) - 1;
          case '0', ordN(ind) = 0;
          case '2', ordN(ind) = max(ordN) + 1;
        end;
%------------------------------------------------------      
%       Reorder eigenvalues using Givens rotations
%------------------------------------------------------      
        while 1
           OK = 1;
           for i=1:n-1
              if ordN(i) > ordN(i+1)
                 if onlyA
                    a11 = A1(i,i); a12 = A1(i,i+1); a22 = A1(i+1,i+1);
                    G1 = cgivens( - a12, a11 - a22 );
                    A1(i:i+1,:) = G1 * A1(i:i+1,:);
                    Q(i:i+1,:)  = G1 * Q(i:i+1,:);
                    A1(:,i:i+1) = A1(:,i:i+1) * G1';
                 else 
                    a11 = A1(i,i); a12 = A1(i,i+1); a22 = A1(i+1,i+1);
                    b11 = B1(i,i); b12 = B1(i,i+1); b22 = B1(i+1,i+1);
   
                    G1 = cgivens( a11*b12-a12*b11, a11*b22-a22*b11 );
                    A1(i:i+1,:) = G1 * A1(i:i+1,:);
                    B1(i:i+1,:) = G1 * B1(i:i+1,:);
                    Q(i:i+1,:)  = G1 * Q(i:i+1,:);
   
                    G2 = cgivens( a12*b22-a22*b12, a11*b22-a22*b11 );
                    A1(:,i:i+1) = A1(:,i:i+1) * G2;
                    B1(:,i:i+1) = B1(:,i:i+1) * G2;
                    Z(:,i:i+1)  = Z(:,i:i+1) * G2;
                 end;
                        
                 t = D(i); D(i) = D(i+1); D(i+1) = t;
                 t = ordN(i); ordN(i) = ordN(i+1); ordN(i+1) = t;
                 OK = 0;
              end
           end
           if OK, break; end;
        end
%------------------------------------------------------      
%       Construct real (generalized) Schur form
%------------------------------------------------------      
        if onlyA
           [Q,A1,D] = csf2rsf (Q', A1);
           B1 = Q;
           Q = D;
        else
           [Q,Z,A1,B1,D] = csf2rsf (Q', Z, triang(A1), triang(B1));
           Q = Q';
           qzcheck(A,A1,Q,Z);
           qzcheck(B,B1,Q,Z);
        end;
         
%------- End of QZPENCIL.M --------- KYuP ----------           
