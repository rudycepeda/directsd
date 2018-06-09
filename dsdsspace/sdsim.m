function [t,X,Y] = sdsim ( sys, K, Tmax )
%SDSIM Impulse response of the standard sampled-data system.
%
%     [T,X,Y] = SDSIM ( SYS, K, TMAX )
%
%   Inputs:
%     SYS - continuous-time generalized plant
%               | z | = | P11 P12 | | w |
%               | y |   | P21 P22 | | u |
%     K   - discrete controller (positive feedback)
%  	  TMAX - final simulation time
%
%   Outputs:
%	  T - time vector
%	  X - state vector
%	  Y - output vector
%
%   See also SDNORM.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%	Check validity
%------------------------------------------------------
        [i2,o2] = size(K);
        [nout,nin] = size(sys);   
        i1 = nin - i2;
        o1 = nout - o2;
        if i1 < 1, error('Negative or zero number of reference inputs'); end;
        if o1 < 1, error('Negative or zero number of working outputs');	end;       
%------------------------------------------------------
%	Prepare to step response simulation instead of
%   delta-impulse: multiply by 's'
%------------------------------------------------------
        sys(:,1:i1) = zpk(sys(:,1:i1)) / tf(1, [1 0]);
        sys = minreal(ss(sys), sqrt(eps), 0);
        n = size(sys.a,1);
        [b1,b2,c1,c2,d11,d12,d21,d22] = blocks4 ( sys.b, sys.c, sys.d, o2, i2 );	
        K = ss(K);
%------------------------------------------------------
%	Simulation using current workspace
%------------------------------------------------------
        options = simset('SrcWorkspace', 'current');
        ws = warning;
        warning off;
        [t,X,Y] = sim('sim_sdstd', Tmax, options);
        warning(ws);
        
%------- End of SDSIM.M --------- KYuP ----------
        
