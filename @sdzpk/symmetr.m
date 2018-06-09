function F = symmetr ( F, type, tol )
%SYMMETR Symmetrize a Hermitian self-conjugate function.
%
%     S = SYMMETR ( F, TYPE, TOL )
%
%   Inputs:
%     F  - initial function
%     TYPE - factorization type of 's' (default), 'z' and 'd'
%     TOL  - relative tolerance (default 1e-2)
%
%   Outputs:
%     S - symmetric fraction
%

%------------------------------------------------------
% Copyright 1999-2006 K. Polyakov 
% $Revision: 3.0$    $Date: 07-Apr-2006 $
%------------------------------------------------------
%       Check data
%------------------------------------------------------
        if ~exist('type','var'), 
          if isct(F), type = 's'; 
          else        type = 'z'; 
          end;
        end;
        if isempty(find(type == 'szd'))
          error('Unknown factorization type ''%s''', type);  
        end;
        if ~exist('tol','var'), tol = 1e-2; end;
%------------------------------------------------------
%       Handle zeros and poles
%------------------------------------------------------
        [z,p,k] = zpkdata ( F, 'v' );
        [zs,zr,z0] = extrpair ( z, type, tol );
        [ps,pr,p0] = extrpair ( p, type, tol );
        g = min(z0,p0);
        z0 = z0 - g; 
        p0 = p0 - g; 
%------------------------------------------------------
%       Brute force symmetrization for large
%       parasitic zeros
%------------------------------------------------------
        if ~isempty(zr)
          ind = find(abs(zr) > 1e3);
          if ~isempty(ind)
            k = k*real(prod(-zr(ind)));
            zr(ind) = [];
          end;
        end;
%------------------------------------------------------
%       Finalization
%------------------------------------------------------
        errmsg = 'Cannot symmetrize the fraction';
        if ~isempty(zr) || ~isempty(pr), error(errmsg); end;
        if type == 's'
          if mod(z0+p0,2) == 1, error(errmsg); end;
          zs = [zs; zeros(z0/2,1)];
          ps = [ps; zeros(p0/2,1)];
          zs = [zs; -zs];
          ps = [ps; -ps];
          F = zpk(zs, ps, k, get(F,'Ts'));           
        else
          pe = length(ps) + p0 - length(zs) - z0;  
          if pe ~= 0, error(errmsg); end;
          zs = [zs; 1./zs; zeros(z0,1)];
          ps = [ps; 1./ps; zeros(p0,1)];
          F = zpk(zs, ps, k, get(F,'Ts'));
        end;

%------- End of SYMMETR.M --------- KYuP ----------
