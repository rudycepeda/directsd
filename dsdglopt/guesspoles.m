function p = guesspoles ( poles, nPoles )
%GUESSPOLES Select N poles as an initial guess.
%
%     PG = GUESSPOLES ( P, N )
%
%   Inputs:
%     P - array of optimal poles
%     N - number of rrequired poles
%
%   Outputs:
%     PG - selected poles with maximal absolute values
%
%   See also MODSDH2, MODSDL2.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Sort by absolute value (decreasing)
%------------------------------------------------------
        [xx,ind] = sort(abs(poles));  
        poles = poles(flipud(ind));
%------------------------------------------------------
%       Select poles
%------------------------------------------------------
        p = zeros(nPoles,1); i = 0;
        while i < nPoles  &&  length(poles) > 0
         %---------------------------------------------
         % Complex pair
         %---------------------------------------------
          if imag(poles(1))~= 0
            if i < nPoles-1  
              i = i + 2;
              p(i-1) = poles(1);
              p(i)   = poles(2);
            end;
            poles(1:2) = [];
         %---------------------------------------------
         % Real pole
         %---------------------------------------------
          else 
            i = i + 1;
            p(i) = poles(1);
            poles(1) = [];
          end;
        end;
%------------------------------------------------------
%       Additional deadbeat poles
%------------------------------------------------------
        p = [p; zeros(nPoles-i,1)];
          
%------- End of GUESSPOLES.M --------- KYuP ----------           
