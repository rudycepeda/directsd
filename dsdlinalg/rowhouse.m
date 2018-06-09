function A = rowhouse ( A, v )
%ROWHOUSE Left (row) Householder transformation.
%
%     B = ROWHOUSE ( A, V )
%
%   Inputs:
%     A - matrix
%     V - vector for Householder transformation
%
%   Outputs:
%     B - transformed matrix B = (I - 2*V*V'/V'*V)*A
%
%   See also HOUSE, COLHOUSE.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        beta = - 2 / (v'*v);
        w = beta * A' * v;
        A = A + v*w';

%------- End of ROWHOUSE.M --------- KYuP ----------           
