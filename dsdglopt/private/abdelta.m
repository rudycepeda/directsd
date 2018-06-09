function [aDelta,bDelta,data] = abdelta ( data, Delta )
%ABDELTA Find polynomials 'aDelta' and 'bDelta'.
%
%	  [ADELTA,BDELTA,DATA] = ABDELTA ( DATA, DELTA )
%
%   Inputs:
%     DATA - structure with the following fields
%         .plant - discrete model of the plant, in 'z' 
%         .dK0   - mandatory part of controller denominator, 
%                  for example [1 -1] to include integrator   
%     DELTA - characteristic polynomial  
%
%   Outputs:
%     ADELTA, BDELTA - minimal controller associated with Delta 
%     DATA - changed DATA structure
%
%   See also GO_PAR2K.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 26-Dec-2006 $
%------------------------------------------------------
%       Construct plant polynomials
%------------------------------------------------------
    if ~isfield(data,'n') ||  ~isfield(data,'d')
       [n,d] = tfdata( data.plant', 'v');
       data.n = striplz(n); data.d = striplz(d);
    end;  
%------------------------------------------------------
%       Find polynomials
%------------------------------------------------------
    dd = conv(data.d,data.dK0);
    if deg(data.n) <= deg(data.d)
         [aDelta,bDelta] = dioph(data.n, dd, Delta);            
    else [bDelta,aDelta] = dioph(dd, data.n, Delta); end;

%------- End of ABDELTA.M --------- KYuP ----------
        