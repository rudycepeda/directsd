function [AB,BB,CB,DB] = bilinss ( sys, K )
%BILINSS State-space bilinear transformation.
%
% 	  SYSB = BILINSS ( SYS, K ) 
% 	  [AB,BB,CB,DB] = BILINSS ( SYS, K ) 
%
%   Inputs:
%	  SYS  - initial system
%	  K    - coefficient 'a' (default K = 1) or a vector of 
%            coefficients [a b c d];
%   Outputs:	
%	  SYSB - bilinear transform of the system SYS, or
%	  AB, BB, CB, DB - its state-space description
%
%   Direct transformation (continuous-time models): 
%            az + b                            - ds + b
%      s = ---------        or            z = ----------  
%            cz + d                             cs - a
%   Inverse transformation (discrete-time models): 
%          - ds + b
%      z = ----------  
%            cs - a
%
%   See also BILINTR.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Find transformation coefficients
%------------------------------------------------------
    	if ~exist('K','var'), K = [1 -1 1 1]; end;
        if length(K) == 1,    K = [K -K 1 1]; end;
    	if length(K) ~= 1  &&  length(K) ~= 4
    	   error('Incorrect coefficient vector');
    	end;
    	if isct(sys), 
    	   a = K(1); b = K(2); c = K(3); d = K(4);
    	   Ts = 1;
        else
    	   a = - K(4); b = K(2); c = K(3); d = - K(1);
    	   Ts = 0;
    	end;
%------------------------------------------------------
%       Conversion
%------------------------------------------------------
        syss = ss(sys); 
        [A,B,C,D] = ssdata(syss);
        n = size(A,1);
        temp = inv(a*eye(n) - c*A);	
        AB = (d*A - b*eye(n)) * temp;
        BB = (a*d - c*b) * temp * B;
        CB = C * temp;
        DB = D + c*C*temp*B;
%------------------------------------------------------
%       Conversion to initial form
%------------------------------------------------------
        if nargout < 2
          sysb = ss ( AB, BB, CB, DB, Ts );    
          if isa(sys, 'tf'),  sysb = tf(sysb); end;
	      if isa(sys, 'zpk'), sysb = zpk(sysb); end;
	      AB = sysb;
        end;

%------- End of BILINSS.M --------- KYuP ----------
