function [bint,bfrac] = val2bin ( x, nfrac )
%VAL2BIN Map a real value to a binary one.
%
%	  [BINT,BFRAC] = VAL2BIN ( X, NF )
%
%   Inputs:
%	  X  - scalar real value
%     NF - bits in fractional part (default 0)
%
%   Outputs:
%	  BINT  - integer part in the form of binary string
%     BFRAC - fractional part in the form of binary string
%
%   See also BIN2VAL.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Integer part
%------------------------------------------------------
        x = abs(x);
        ix = floor(x);
        fx = x - ix;
        bint = dec2bin ( ix );
%------------------------------------------------------
%       Fractional part
%------------------------------------------------------
        if ~exist('nfrac','var'), nfrac = 0; end;
        bfrac = [];
        for i=1: nfrac
          fx = fx * 2;
          ix = floor ( fx );
          fx = fx - ix;
          bfrac = [bfrac char(48+ix)];
        end;
        
%------- End of VAL2BIN.M --------- KYuP ----------
