function v = house ( x, xx )
%HOUSE Householder vector V for a given vector X.
%
%     V = HOUSE ( X )
%     V = HOUSE ( X, 1 )
%
%   Inputs:
%     X - column-vector
%
%   Outputs:
%     V - vector such that 
%           (I-2*v*v'/v'*v) x = [* 0 ... 0]'
%         or, with specified second parameter,
%           x' (I-2*v*v'/v'*v) = [0 ... 0 *]
%
%   See also COLHOUSE, ROWHOUSE.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Preliminaries
%------------------------------------------------------
        n = length(x);
        mu = norm(x);
        v = x;
%------------------------------------------------------
%       Normal case
%------------------------------------------------------
        if exist('xx','var')
           if mu ~= 0
              if x(n) < 0, signum = -1;
              else         signum = 1; end;
              beta = x(n) + signum*mu;
              v(1:n-1) = v(1:n-1) / beta;
           end;
           v(n) = 1;
%------------------------------------------------------
%       Special case
%------------------------------------------------------
        else
           if mu ~= 0
              if x(1) < 0, signum = -1;
              else         signum = 1; end;
              beta = x(1) + signum*mu;
              v(2:n) = v(2:n) / beta;
           end;
           v(1) = 1;
        end;

%------- End of HOUSE.M --------- KYuP ----------           
