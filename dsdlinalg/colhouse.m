function A = colhouse ( A, v )
%COLHOUSE Right (column) Householder transformation.
%
%     B = COLHOUSE ( A, V )
%
%   Inputs:
%     A - matrix
%     V - vector for Householder transformation
%
%   Outputs:
%     B - transformed matrix B = A*(I - 2*V*V'/V'*V)
%
%   See also HOUSE, ROWHOUSE.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        beta = - 2 / (v'*v);
        w = beta * A * v;
        A = A + w*v';
        
%------- End of COLHOUSE.M --------- KYuP ----------           
