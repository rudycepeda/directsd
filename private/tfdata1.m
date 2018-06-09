function [n,d] = tfdata1 ( F )
%TFDATA1 Extract numerator and denominator removing leading zeros.
%
% 	  [N,D] = TFDATA1 ( F )
%
%   Inputs:
%	  F - an LTI model
%
%   Outputs:
%	  N, D - polynomial numerator and denominator
%
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
	    [n,d] = tfdata ( F, 'v' );
	    n = striplz ( n );
	    d = striplz ( d );

%------- End of TFDATA1.M --------- KYuP ----------
