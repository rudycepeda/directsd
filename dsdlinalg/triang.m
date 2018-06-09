function A = triang ( A )
%TRIANG Zeroing all elements below main diagonal.
%
%     B = TRIANG ( A )
%
%   Inputs:
%     A - matrix supposed to be in the upper 
%         triangular form
%
%   Outputs:
%     B - matrix in uper triangular form
%
%   See also HESSFORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        n = size(A, 1);
        tol = 100*eps;
%------------------------------------------------------
        for i=1:n
           for j=1:i-1
              A(i,j) = 0;
           end;
           for j=i:n
              if abs(A(i,j)) < tol, A(i,j) = 0; end;              
              if abs(imag(A(i,j))) < tol,
                 A(i,j) = real(A(i,j));
              end;              
           end;           
        end;

%------- End of TRIANG.M --------- KYuP ----------           
