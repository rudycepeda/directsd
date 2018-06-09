function [rho] = cp2par ( DeltaZ, alpha, beta, nPairs )
%CP2PAR Map a polynomial to an array of parameters in [0,1]. 
%
%	  RHO = CP2PAR ( DZ, ALPHA, BETA, NPAIRS )
%
%   Inputs:
%     DZ     - characteristic polynomial in z
%     ALPHA  - limiting stability degree
%     BETA   - limiting damping factor
%                BETA>0 for truncated sector and
%                BETA<0 for shifted sector
%     NPAIRS - number of complex conjugated pairs (other poles are real)  
%              (default - maximal possible)
%   
%   Outputs:
%     RHO    - vector of parameters in [0,1] 
%
%   See also PAR2CP.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('nPairs','var'), nPairs = floor(deg(DeltaZ)/2); end;
        if ~exist('alpha','var'),  alpha = 0; end;
        if ~exist('beta','var'),   beta = Inf; end;
%------------------------------------------------------
%       Find roots of DeltaZ and project them
%       onto the admissible region
%------------------------------------------------------
        poles = roots(DeltaZ);
        poles = admproj(poles, alpha, beta, 'z');
        [sp,ind] = sort(abs(imag(poles)));
        poles = poles(flipud(ind));
        nP = length(find(sp > 0)) / 2;
        if nP > nPairs,
          error('The polynomial has more complex pairs than expected'); 
        end;
%------------------------------------------------------
%       Check syntax: beta<0 for shifted sector
%------------------------------------------------------
        Shifted = (beta < 0);
        beta = abs(beta);
        Ea = exp(-alpha);
        Eb = exp(-pi/beta);
        if ~Shifted  &&  Eb > Ea
          beta = Inf; Eb = 1;  
        end;
%------------------------------------------------------
%       Limits for r2
%------------------------------------------------------
        [r2min,r2max,E0] = r2range(Ea, Eb, Shifted);
%------------------------------------------------------
%       Paired
%------------------------------------------------------
        rho = poles*0;
        k = 1;
        for i=1:nPairs
          p = poly(poles(k:k+1));  
          r1 = p(2); r2 = p(3);
          rho(k) = (r2 - r2min)/(r2max - r2min);
          [r1min,r1max] = r1range ( r2, Ea, beta, Shifted );
          rho(k+1) = (r1 - r1min)/(r1max - r1min);
          k = k + 2;
        end;
%------------------------------------------------------
%       First-degree terms
%------------------------------------------------------
        r0min = - Ea;
        r0max = E0; 
        while k<=length(poles)
          r0 = - poles(k);  
          rho(k) = (r0 - r0min)/(r0max - r0min);
          k = k + 1;
        end;

%------- End of CP2PAR.M --------- KYuP ----------
