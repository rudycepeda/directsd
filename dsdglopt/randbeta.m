function rb = randbeta ( a, b, rows, cols )
%RANDBETA Beta-distributed random numbers.
%
%       X = RANDBETA ( A, B, R, C )
%
%   Inputs: 
%       A, B - parameters
%       R - rows
%       C - columns
%
%   Outputs:
%       X - array of beta-distributed values
%
%   See also RANDGAMMA.

% [1]  Devroye, L. Non-Uniform Random Variate Generation, 
%      Springer-Verlag, 1986.
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%%       Check data
%------------------------------------------------------
        if ~exist('rows','var'), rows = 1; end;
        if ~exist('cols','var'), cols = rows; end;
%------------------------------------------------------
%%      Generate gamma random values and take ratio 
%%      of the first to the sum.
%------------------------------------------------------
        g1 = randgamma(a, rows, cols);
        g2 = randgamma(b, rows, cols);
        rb = g1 ./ (g1 + g2);

%------- End of RANDBETA.M --------- KYuP ----------
