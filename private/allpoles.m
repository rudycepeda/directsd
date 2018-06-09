function P = allpoles ( sys )
%ALLPOLES Find all poles of zpk-matrix.
%       
%     P = ALLPOLES ( SYS )
%
%   Inputs:
%     SYS - a square LTI system model in zpk-form 
%
%   Outputs:
%     P - all poles counting multiplicities
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        [rows,cols] = size(sys);        
        P = [];        
        S.type = '()';
        for r=1:rows
          for c=1:cols
             S.subs = {r, c};  
             F = subsref(sys, S);
             p = others(F.p{1}, P, sqrt(eps));
             P = [P; p];
          end;
        end;
        
%------- End of ALLPOLES.M --------- KYuP ----------           
