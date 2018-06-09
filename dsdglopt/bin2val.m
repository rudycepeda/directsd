function x = bin2val ( bint, bfrac  )
%BIN2VAL Transform a binary string value to a real number.
%
%	  X = BIN2VAL ( BINT, BFRAC )
%
%   Inputs:
%	  BINT  - integer part in the form of binary string
%     BFRAC - fractional part in the form of binary string
%
%   Outputs:
%	  X  - scalar real value
%
%   See also VAL2BIN.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Integer part
%------------------------------------------------------
        x = bin2dec ( bint );
        lfrac = length(bfrac);
        x = x + bin2dec(bfrac) / 2^lfrac;
%------- End of VAL2BIN.M --------- KYuP ----------

