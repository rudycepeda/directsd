function [r1min,r1max,E0] = r1range ( r2, Ea, beta, Shifted )
%R1RANGE Find the range of r1.
%
%	  [R1MIN,R1MAX,E0] = R1RANGE ( R2, EA, BETA, SHIFTED )
%
%   Inputs:
%     R2 - value of free term r2 
%     EA = exp(-\alpha*T), where \alpha denotes degree of stability
%     BETA denotes oscillation
%     SHIFETD = 0 for a truncated sector
%             = 1 for a shifted sector  
%   
%   Outputs:
%     R1MIN, R1MAX - lower and upper limits for r2 
%     E0 - auxiliary value
%
%   See also R2RANGE.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        Eb = exp(-pi/beta);
%------------------------------------------------------
        if ~Shifted           % truncated sector   
             E0 = min(Ea,Eb);    
        else E0 = Ea*Eb;      % shifted sector 
        end; 
%------------------------------------------------------
        r2crit = E0^2;
        r1min = -r2/Ea - Ea;
        if r2 <= r2crit
           r1max = r2/E0 + E0; 
        else 
           sqr2 = sqrt(r2); 
           if Shifted, sqr2 = sqr2 / Ea; end;
           r1max = -2*sqrt(r2)*cos(-beta*log(sqr2)); 
        end;

%------- End of R1RANGE.M --------- KYuP ----------
