function [Delta,D22] = charpol ( sys, K, H )
%CHARPOL Characteristic polynomial of sampled-data system.
%
%     P = CHARPOL ( SYS, K, H )
%
%   Inputs:
%     SYS - a standard LTI plant
%     K   - an LTI  controller
%     H   - an extrapolator (default is ZOH)
%
%   Outputs:
%     P   - characteristic polynomial of the closed loop 
%
%   See also SDMARGIN.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Construct coefficients
%------------------------------------------------------
        sys = zpk(sys);
        if K.Ts ~= sys.Ts,
           P22 = sys(end,end);
           if nargin <= 2
                D22 = dtfm(P22, K.Ts, 0);
           else D22 = dtfm(P22, K.Ts, 0, H); 
           end;
        else
           D22 = sys(end,end);
        end;
        D22 = - D22;
        X = 1 + D22*K;
        Delta = tf2nd ( X );

%------- End of CHARPOL.M --------- KYuP ----------           
        
