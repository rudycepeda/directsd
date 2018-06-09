function x = u2range ( u, min, max )
%U2RANGE Linear mapping from [0,1] onto [min,max].
%
%	  X = U2RANGE ( U, MIN, MAX )
%
%   Inputs:
%	  U - a value such that 0 <= U <= 1 
%	  MIN, MAX - range bounds
%
%   Outputs:
%	  X - image
%
%   See also UNIPROJ.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        x = min + (max - min)*u;

%------- End of U2RANGE.M --------- KYuP ----------
