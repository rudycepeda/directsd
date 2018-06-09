function [Q,A,R] = compress ( A, type, tol )
%COMPRESS Row or column compression for a matrix.
%         Method: left (or right) multiplication by a 
%         unitary matrix Q such that Q*Q'=I.
%
%   Row compression:
%  ------------------
%     [Q,B,R] = COMPRESS ( A, 'row' )
%     [Q,B,R] = COMPRESS ( A, 'row', TOL )
%
%   Inputs:
%     A   - initial matrix
%     TOL - tolerance for rank evaluation
%
%   Outputs:
%     Q - unitary matrix such that Q*Q'=Q'*Q=I
%     B - a compressed matrix with zero bottom rows
%         such that Q * A = B  
%     R - rank of B and the number of nonzero rows 
%
%   Column compression:
%  ---------------------
%     [Z,B,R] = COMPRESS ( A, 'col' )
%     [Z,B,R] = COMPRESS ( A, 'col', TOL )
%
%   Inputs:
%     A   - initial matrix
%     TOL - tolerance for rank evaluation
%
%   Outputs:
%     Z - unitary matrix such that Z*Z'=Z'*Z=I
%     B - a compressed matrix with zero left columns
%         such that A * Z = B
%     R - rank of B and the number of nonzero columns
%
%   See also QR.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('tol', 'var'), tol = 1.e-6; end;       
        if ~isequal(type,'row') && ~isequal(type,'col')
           error('Unknown transformation type');
        end;
        if isequal(type,'col'), A = A'; end;
%------------------------------------------------------
%       Rank evaluation using QR-decomposition
%------------------------------------------------------
        [n,m] = size(A);       
        [Q,A,E] = qr(A);
        A = A*E';
        R = rank(A, tol);
        if R < n
           A(R+1:n,:) = zeros(n-R,m);
        end;
%------------------------------------------------------
%       Finalization
%------------------------------------------------------
        if isequal(type,'col')
          A = A';
          Q = Q';
          if R < n
             ind = [R+1:n 1:R];
             A = A(:,ind);
             Q = Q(ind,:);
          end;
        end;       
        Q = Q';

%------- End of COMPRESS.M --------- KYuP ----------           
