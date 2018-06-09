function K = optkdelta ( data, Delta, degKsi, gKsi, aDelta, bDelta )
%OPTKDELTA Find optimal Delta-controller.
%
%	  K = OPTKDELTA ( DATA, DELTA, DEG, GAIN )
%
%   Inputs:
%     DATA - structure with the following fields
%         .A, .B - coefficients in the functional, in 'zeta'
%         .plant - discrete model of the plant, in 'z' 
%         .dK0   - mandatory part of controller denominator, 
%                  for example [1 -1] to include integrator   
%     DELTA - characteristic polynomial  
%     DEG   - degree of the 'ksi'-polynomial  
%     GAIN  - 'static gain' of the 'ksi'-polynomial (NaN if any) 
%
%   Outputs:
%     K - optimal controller associated with Delta 
%
%   See also OPTKSI, GO_PAR2K.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 26-Dec-2006 $
%------------------------------------------------------
%       Check params
%------------------------------------------------------
    if ~exist('gKsi','var'), gKsi = NaN; end;
    if ~isfield(data,'dK0'), data.dK0 = 1; end;
%------------------------------------------------------
%       Generate aDelta & bDelta
%------------------------------------------------------
    if ~exist('aDelta','var') || ~exist('bDelta','var')
      [aDelta,bDelta,data] = abdelta ( data, Delta );  
    end;
%------------------------------------------------------
%       Find optimal ksi
%------------------------------------------------------
    ksi = optksi ( data, Delta, degKsi, gKsi, aDelta, bDelta );
%------------------------------------------------------
%       Construct controller
%------------------------------------------------------
    num = sumpol(aDelta,  conv(conv(data.d,ksi), data.dK0));
    den = conv(sumpol(bDelta, -conv(data.n,ksi)), data.dK0);
    K = zpk(tf(num, den, data.plant.Ts)');
        
%------- End of OPTKDELTA.M --------- KYuP ----------
