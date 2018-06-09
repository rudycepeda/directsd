function c = sumpol (a, b, varargin )
%SUMPOL Sum of two polynomials of different degrees.
%
%     C = SUMPOL ( A, B )
%     
%   Inputs:
%     A, B - polynomials as coefficient arrays
%
%   Outputs:
%     C - the sum as an array
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Multiple polynomials version
%------------------------------------------------------
        if nargin > 2, 
           c = sumpol(a, b);
           for i=3:nargin
             c = sumpol(c, varargin{i-2});  
           end;
           return;
        end;
%------------------------------------------------------
%       Simplest version for two polynomials
%------------------------------------------------------
        la = length(a);
        lb = length(b);
        c = [zeros(1,lb-la) a] + [zeros(1,la-lb) b];

%------- End of SUMPOL.M --------- KYuP ----------
