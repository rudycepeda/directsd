function [n,d] = tf2nd ( F, var )
%TF2ND Extract numerator and denominator polynomials.
%
%     [N,D] = TF2ND ( F, VAR )
%
%   Inputs:
%     F   - LTI model (in 'ss', 'tf' or 'zpk' form)
%     VAR - variable for polynomials (default is 's' for  
%           continuous-time models and 'z' for discrete-time ones
%       
%   Outputs:
%     N, D - numerator and denominator polynomials
%
%   See also IMPROPER.

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.00 $    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Determine variable
%------------------------------------------------------
        if ~exist('var','var'), 
          if isct(F)
              var = 's'; 
          else
              var = 'z';
          end
        end
%------------------------------------------------------
%       Extraction
%------------------------------------------------------
        if isa(F,'zpk')
          [z,p,k] = zpkdata(F,'v');
          n = k*poln(z, ['r' var]);
          d = poln(p, ['r' var]);
        else
          [n,d] = tfdata(F, 'v');
          n = poln(n, var);
          d = poln(d, var);
        end;

%------- End of TF2ND.M --------- KYuP ----------
