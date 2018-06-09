function g = sdmargin ( sys, K, H, type )
%SDMARGIN Stability margin for sampled-data system.
%
%     G = SDMARGIN ( SYS, K, H, TYPE )
%
%   Inputs:
%     SYS  - an LTI system
%     K    - discrete-time controller
%     H(s) - transfer matrix of the hold:
%              h(t)=K*exp(L*t)*M => H(s)=K*(sI-L)^{-1}*M
%            default hold is the scalar zero-order hold (ZOH) with H(s)=1/s.
%     TYPE - type of margin: 'hinf' of 'ahinf' (default)
%
%   Outputs:
%     G  - stability margin with respect to periodically 
%	       nonstationary perturbations	
%
%   See also SDAHINORM, SDHINORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
        if ~exist('H','var'), H = []; 
        elseif ischar(H)
          type = H; H = [];  
        end;
        if ~exist('type','var') || isempty(type) 
           type = 'ahinf'; 
        end;
        switch type 
          case 'hinf',  g = 1 / sdhinorm(sys, K);
          case 'ahinf', g = 1 / sdahinorm(sys, K, H);
        otherwise
          error('Unknown type of stability margin ''%s''', type);  
        end;

%------- End of SDMARGIN.M --------- KYuP ----------           
