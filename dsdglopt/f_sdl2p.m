function [f,K] = f_sdl2p ( coef )
%F_SDL2P Target function for nonlinear L2-optimization.
%
%	  [F,K] = F_SDL2P ( COEF )
%
%   Inputs:
%     COEF - vector of parameters
%   
%   Outputs:
%	  F - value of L2-error
%     K - controller associated with X
%
%   See also F_SDH2P.
	
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
%       Compute integral error
%------------------------------------------------------
        %f = sdl2err (gloptData.sys, K);        
       %-----------------------------------------------
       %  Use lifted model
       %-----------------------------------------------
        if gloptData.ss            
          ws = warning; warning off;  
          Fcl = minreal(lft(gloptData.dsys,K),sqrt(eps),0);  
          f = K.Ts*norm(Fcl)^2;
          if isinf(f)  ||  f == 0
            Fcl = minreal(zpk(lft(gloptData.dsys,K)),1e-3);  
            f = K.Ts*norm(Fcl)^2;
            if isinf(f)
              %disp('Infinite value of the L2-cost function');  
              %keyboard;
            end;
          end;
          warning(ws);
       %-----------------------------------------------
       %  Use frequency-domain formulas
       %-----------------------------------------------
        else
          M = feedback ( -K, -gloptData.plant )';
          f = quaderr ( gloptData.A, -gloptData.B, gloptData.E, M ); 
        end;
        
%------- End of F_SDL2P.M --------- KYuP ----------
