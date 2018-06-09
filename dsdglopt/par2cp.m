function [Delta,DeltaZ,poles] = par2cp ( rho, alpha, beta, nPairs )
%PAR2CP Map array of parameters in [0,1] to the characteristic polynomial.
%
%	  [D,DZ] = PAR2CP ( RHO, ALPHA, BETA, NPAIRS )
%
%   Inputs:
%     RHO    - vector of parameters in [0,1] 
%     ALPHA  - limiting stability degree
%     BETA   - limiting damping factor
%              BETA>0 for truncated sector and
%              BETA<0 for shifted sector
%     NPAIRS - number of complex conjugated pairs (other poles are real)  
%              (default - maximal possible)
%   
%   Outputs:
%     D  - characteristic polynomial in zeta
%     DZ - characteristic polynomial in z
%
%   See also CP2PAR.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('nPairs','var'), nPairs = floor(length(rho)/2); end;
        if ~exist('alpha','var'),  alpha = 0; end;
        if ~exist('beta','var'),   beta = Inf; end;
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
        poles = rho*0; k = 1;
        DeltaZ = 1;
        for i=1:nPairs
          r2 = r2min + rho(k)*(r2max - r2min);
          [r1min,r1max] = r1range ( r2, Ea, beta, Shifted );
          r1 = r1min + rho(k+1)*(r1max - r1min);
          d = [1 r1 r2];
          poles(k:k+1) = roots(d);
          DeltaZ = conv(DeltaZ, d);
          k = k + 2;
        end;
%------------------------------------------------------
%       Single
%------------------------------------------------------
        r0min = - Ea;
        r0max = E0; 
        while k<=length(rho)
          r0 = r0min + rho(k)*(r0max - r0min);
          if abs(r0) < 1e-3, r0 = 0; end;
          poles(k) = - r0;
          DeltaZ = conv(DeltaZ, [1 r0]);
          k = k + 1;
        end;
%------------------------------------------------------
%       Reciprocal polynomial
%------------------------------------------------------
        Delta = striplz(fliplr(DeltaZ));

%------- End of PAR2CP.M --------- KYuP ----------
