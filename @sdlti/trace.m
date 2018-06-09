function T = trace ( sys )
%TRACE Trace of a square system model.
%       
%     T = TRACE ( SYS )
%
%   Inputs:
%     SYS - a square LTI system model in zpk-form 
%
%   Outputs:
%     T - trace (the sum of diagonal elements)
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        [rows,cols] = size(sys);        
        if rows ~= cols, 
           error('Trace is defined only for square matrices');
        end;
%------------------------------------------------------
%       Find trace in 'ss' or 'zpk' form
%------------------------------------------------------
        S.type = '()';
        T = tf(0, 1, sys.Ts);
        if isa(sys,'ss')
             T = ss(T);
        else T = zpk(T); 
        end;
        for i=1:rows
          S.subs = {i, i};
          if isa(sys,'ss')
               T = minreal(T + subsref(sys,S), sqrt(eps), 0);
          else T = minreal(sumzpk(T,subsref(sys,S))); 
          end;
        end;
%------------------------------------------------------
%       Final form
%------------------------------------------------------
        if isa(sys,'tf')
           T = tf(T); 
        end;
        
%------- End of TRACE.M --------- KYuP ----------           
