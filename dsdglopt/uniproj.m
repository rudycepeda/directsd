function x = uniproj ( x )
%UNIPROJ Project values onto [0,1].
%
%   See also U2RANGE.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        x = max(0, min(x,1));
   
%------- End of UNIPROJ.M --------- KYuP ----------   
