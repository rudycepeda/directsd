function [r2min,r2max,E0] = r2range ( Ea, Eb, Shifted )
%R2RANGE Find the range of r2.
%
%	  [R2MIN,R2MAX,E0] = R2RANGE ( EA, EB, SHIFTED )
%
%   Inputs:
%     EA = exp(-\alpha*T)
%     EB = exp(-\pi/\beta)
%     SHIFETD = 0 for a truncated sector
%             = 1 for a shifted sector  
%   
%   Outputs:
%     R2MIN, R2MAX - lower and upper limits for r2 
%     E0 - auxiliary value
%
%   See also R1RANGE.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        r2max = Ea^2;
        if ~Shifted   % truncated sector   
          E0 = min(Ea,Eb);    
          r2min = -Ea*E0;
        else          % shifted sector 
          r2min = -Ea^2*Eb;
          E0 = Ea*Eb; 
        end        

%------- End of R2RANGE.M --------- KYuP ----------
