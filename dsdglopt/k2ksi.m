function [ksi,aDelta,bDelta] = k2ksi ( sys, K, dK0, H )
%K2KSI Find a polynomial parameter \xi for a controller.
%
%	  KSI = K2KSI ( SYS, K, DK0 )
%	  KSI = K2KSI ( SYS, K, DK0, H )
%
%   Inputs:
%     SYS - continuous-time plant in standard system
%     K   - discrete-time controller
%     DK0 - mandatory part of controller denominator 
%     H   - generalized hold as an LTI model 
%
%   Outputs:
%     KSI - polynomial parameter such that
%             K = (aDelta + d*ksi)/(bDelta - n*ksi)
%
%   See also PAR2CP, CP2PAR.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        T = K.Ts;
        if ~exist('dK0','var'), dK0 = 1; end;
        if ~exist('H','var'), H = []; end;
%------------------------------------------------------
%       Find characteristic polynomial
%------------------------------------------------------
        [DeltaZ,Pz] = charpol ( sys, K, H );
        Delta = recip(DeltaZ);
        Delta = Delta.coef;
%------------------------------------------------------
%       Find aDelta and bDelta 
%------------------------------------------------------
        [n,d] = tfdata( Pz', 'v');
        n = striplz(n); d = striplz(d);
        dd = conv(d,dK0);
        if deg(n) <= deg(d)
             [aDelta,bDelta] = dioph(n, dd, Delta);            
        else [bDelta,aDelta] = dioph(dd, n, Delta); 
        end;
%------------------------------------------------------
%       Find ksi
%------------------------------------------------------
        K = K';
        ksiN = sumzpk(minreal(K*tf(bDelta,1,T)),-tf(aDelta,1,T));
        ksiD = sumzpk(minreal(K*tf(n,1,T)),tf(d,1,T));
        ksi = minreal(ksiN/ksiD);
        [ksi,dKsi] = tf2nd(ksi);
        [ksi,err] = mrdivide(ksi,dKsi);
        if norm(err) > 1e-3*norm(dKsi),
          disp(sprintf('Cancellation error %g',norm(err)));  
        end;

%------- End of K2KSI.M --------- KYuP ----------
