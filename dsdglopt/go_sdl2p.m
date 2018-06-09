function [f,coef,K] = go_sdl2p ( x )
%GO_SDL2P Target function for modal-based L2-optimization.
%
%	  [F,COEF,K] = GO_SDL2P ( X )
%
%   Inputs:
%	  X - scalar optimization parameter in [0,1] 
%   
%   Outputs:
%	  F - constraints values and function value to be minimized 
%     K - controller associated with X
%
%   See also GO_SDH2P.
	
%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
        global gloptData
        if isobject(x) % Special case: x = controller
             coef = x;
        else coef = hilb2coord ( x, gloptData.dim, gloptData.bits );
        end;
%------------------------------------------------------
%       Find cost function
%------------------------------------------------------
        [f,K] = f_sdl2p ( coef );
        
%------- End of GO_SDL2P.M --------- KYuP ----------
