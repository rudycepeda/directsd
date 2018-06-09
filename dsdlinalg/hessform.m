function A = hessform ( A )
%HESSFORM Zeroing all elements of a matrix below the first subdiagonal.
%         In the Hessenberg form they must be zero.
%
%     B = HESSFORM ( A )
%
%   Inputs:
%     A - matrix supposed to be in the Hessenberg form
%
%   Outputs:
%     B - matrix in Hessenberg form
%
%   See also TRIANG.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
     n = size(A, 1);
     for i=3:n, for j=1:i-2, A(i,j)=0; end; end;
     for i=2:n
        if abs(A(i,i-1)) < 1.e-10*(abs(A(i-1,i-1)) + abs(A(i,i)))
           A(i,i-1) = 0;
        end;
     end;

%------- End of HESSFORM.M --------- KYuP ----------           
