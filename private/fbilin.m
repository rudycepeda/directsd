function BL = fbilin ( b )
%FBILIN Bilinear transformation factor.
%
%     F = BILIN ( B )
%
%   Inputs:
%     B - a constant (default 1)
%
%   Outputs:
%     F - additional factor due to non-standard 
%         bilinear transform
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Additional term in the integral due to
%       nonstandard bilinear transformation
%              z = (s + 1) / (s - 1)
%------------------------------------------------------
%                     - 2 b 
%	Hence, dz/z = -------------- ds
%                 (s + b)(s - b)
%
%------------------------------------------------------
        if ~exist('b','var'), b = 1; end;
        BL = zpk([], [b -b], -2*b);
        
%------- End of BILIN_0.M --------- KYuP ----------
