function [f,K] = f_sdh2p ( coef )
%F_SDH2P Target function for nonlinear H2-optimization.
%
%	  [F,K] = F_SDH2P ( COEF )
%
%   Inputs:
%     COEF - vector of parameters
%   
%   Outputs:
%	  F - value of H2-norm
%     K - controller associated with X
%
%   See also F_SDL2P.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        global gloptData
        if isobject(coef)
             K = coef; 
        else K = go_par2k ( coef ); 
        end;
%------------------------------------------------------
%       Compute H2-norm
%------------------------------------------------------
        %f = sdh2norm (gloptData.sys, K);        
       %-----------------------------------------------
       %  Use lifted model
       %-----------------------------------------------
        if gloptData.ss
          ws = warning; warning off;  
          Fcl = minreal(lft(gloptData.dsys,K),sqrt(eps),0);  
          err = norm(Fcl)^2;
          if isinf(err) || err == 0
            Fcl = minreal(zpk(lft(gloptData.dsys,K)),1e-3);  
            err = norm(Fcl)^2;
            if isinf(err)
              %disp('Infinite value of the L2-cost function');  
              %keyboard;
            end;
          end;
          warning(ws);
          err = err + gloptData.gamma;  
       %-----------------------------------------------
       %  Use frequency-domain formulas
       %-----------------------------------------------
        else
           M = feedback ( -K, -gloptData.plant )';
           err = quaderr ( gloptData.A, -gloptData.B, gloptData.E, M ); 
        end;
        f = sqrt(err);
        
%------- End of F_SDH2P.M --------- KYuP ----------
