function K = go_par2k ( coef )
%GO_PAR2K Map array of parameters in [0,1] to the controller.
%
%	  K = GO_PAR2K ( RHO )
%
%   Inputs:
%     RHO - vector of parameters in [0,1]
%
%   Outputs:
%     K   - controller associated with RHO
%
%   See also PAR2CP, OPTKSI, OPTKDELTA.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 26-Dec-2006 $
%------------------------------------------------------
%       Construct characteristic polynomial
%------------------------------------------------------
        global gloptData
        Delta = par2cp ( coef, gloptData.alpha, ...
                         gloptData.beta, gloptData.nPairs );
%------------------------------------------------------
%       Construct aDelta and bDelta
%------------------------------------------------------
        [aDelta,bDelta,gloptData] = abdelta ( gloptData, Delta );  
%------------------------------------------------------
%       Find degree of ksi with acccount for static gain
%------------------------------------------------------
        degKsi = gloptData.dim - 2*deg(gloptData.plant) ...
                 - deg(gloptData.dK0);
        if ~isnan(gloptData.gPsi), 
           degKsi = degKsi - 1; 
           gKsi = polyval(gloptData.a0,1)*polyval(bDelta,1) - ...
                  polyval(gloptData.b0,1)*polyval(aDelta,1) + ...
                  polyval(Delta,1)*gloptData.gPsi; 
        else
           gKsi = NaN;
        end;
        if degKsi < -1,
          disp('Warning: Controller order is too low for such a plant');  
          disp(sprintf('         Order %d is taken.',gloptData.ordK-1-degKsi));
        end;
%------------------------------------------------------
%       Construct controller
%------------------------------------------------------
        K = optkdelta ( gloptData, Delta, degKsi, gKsi, aDelta, bDelta );
        
%------- End of GO_PAR2K.M --------- KYuP ----------
