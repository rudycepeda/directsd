function s = expandstr(s, newLen, b)
%EXPANDSTR Expand string symmetrically.
%
%     SE = EXPANDSTR ( S, LEN, W )
%
%   Inputs:
%     S   - initial string
%     LEN - new length
%     W   - symbol for expansion  
%
%   Outputs:
%     SE - expanded string
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        len = length(s);
        right = fix( (newLen - len) / 2 );
        left = newLen - len - right;
        s = [b(ones(1,left)) s b(ones(1,right))];

%------- End of EXPANDSTR.M --------- KYuP ----------
