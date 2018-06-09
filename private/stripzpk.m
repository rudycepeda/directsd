function F = stripzpk ( F, tol )
%STRIPZPK Eliminate fanthom large zeros of zpk-model.
%	      These very large zeros are caused by summation 
%         of models or similar non-robust operations.
%
%     G = STRIPZPK ( F [, TOL] )
%
%   Inputs:
%	  F   - initial model in zpk-form
%     TOL - tolerance (default 1.e-6)
%
%   Outputs:
%     G   - transformed model
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data    
%------------------------------------------------------
        if ~exist('tol','var'), tol = 1.e-6; end;
        tol = 1 / tol;
%------------------------------------------------------
%       MIMO-case
%------------------------------------------------------
        if max(size(F)) > 1,
          for r=1:size(F,1)  
            for c=1:size(F,2)  
              F(r,c) = stripzpk(F(r,c));  
            end;
          end;
          return;  
        end;
%------------------------------------------------------
%       Remove leading numerator zeros    
%------------------------------------------------------
        n = striplz(F.K*poly(F.z{1}), eps);
        F.z{1} = roots(n);
        F.K = n(1);
%------------------------------------------------------
%       Remove phantom zeros    
%------------------------------------------------------
        tol2 = eps;
        while find(abs(F.z{1}) > tol)
          n = striplz(F.K * poly(F.z{1}), tol2);
          F.z{1} = roots(n);
          F.K = n(1);
          tol2 = tol2 * 10;
        end;
                
%------- End of STRIPZPK.M --------- KYuP ----------
