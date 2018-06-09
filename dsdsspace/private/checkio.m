function  [o1,i1] = checkio ( nout, nin, o2, i2, canBeZero )
%CHECKIO Check input and output dimensions
%
%      [O1,I1] = CHECKIO ( NOUT, NIN, NMEAS, NCON )
%
%   Inputs:
%     NOUT  - number of outputs
%     NIN   - number of inputs
%     NMEAS - dimension of the vector 'y' (default 1)
%     NCON  - dimension of the vector 'u' (default 1)
%
%   Outputs:
%     O1 - number of external outputs
%	  I1 - number of external outputs
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
    if ~exist('canBeZero','var'), canBeZero = 0; end;
	i1 = nin - i2;
	o1 = nout - o2;
    if canBeZero
      if i2 < 0, error('Negative number of control signals'); end;
  	  if o2 < 0, error('Negative number of  measured signals'); end;
    else
  	  if i2 < 1,  error('No control signals specified'); end;
	  if o2 < 1, error('No measured signals specified'); end;       
    end;
	if i1 < 1, error('Negative or zero number of reference inputs'); end;
	if o1 < 1, error('Negative or zero number of working outputs');	end;

%------- End of CHECKIO.M --------- KYuP ----------
